import os
import tempfile
import csv
import json
import collections
import numpy as np

from ..shape import InputShape
from ..shape import InputShapeSingle, InputShapeList, InputShapePairOfLists

from ... import logger

MIN_COLUMN_VALIDITY = 0.8
FLATTENED_EVIDENCE = 0.2
DEFAULT_DELIMITER = ","
FILE_CHUNKSIZE = 10000


class FileTyper(object):
    def __init__(self, path):
        self.path = os.path.join(path)

    def is_valid_input_file(self):
        if self.is_csv() or self.is_tsv() or self.is_json():
            return True
        else:
            return False

    def is_valid_output_file(self):
        if self.is_csv() or self.is_tsv() or self.is_json() or self.is_hdf5():
            return True
        else:
            return False

    def is_tabular(self):
        if self.is_csv() or self.is_tsv():
            return True
        else:
            return False

    def is_csv(self):
        if self.path.endswith(".csv"):
            return True
        else:
            return False

    def is_tsv(self):
        if self.path.endswith(".tsv"):
            return True
        else:
            return False

    def is_hdf5(self):
        if self.path.endswith(".h5"):
            return True
        else:
            return False

    def is_json(self):
        if self.path.endswith(".json"):
            return True
        else:
            return False

    def get_extension(self):
        return self.path.split(".")[-1]


class BatchCacher(object):
    def __init__(self):
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-")

    def get_cached_files(self, prefix):
        idx2fn = {}
        for fn in os.listdir(self.tmp_folder):
            if fn.startswith(prefix):
                idx = int(fn.split("-")[-1].split(".")[0])
                idx2fn[idx] = os.path.join(self.tmp_folder, fn)
        idxs = sorted(idx2fn.keys())
        return [idx2fn[idx] for idx in idxs]

    def get_cached_input_files(self):
        return self.get_cached_files(prefix="chunk-input-")

    def get_cached_output_files(self):
        return self.get_cached_files(prefix="chunk-output-")

    def name_cached_output_files(self, cached_inputs, output_template):
        ft = FileTyper(output_template)
        extension = ft.get_extension()
        cached_outputs = []
        for i, _ in enumerate(cached_inputs):
            cached_outputs += [
                os.path.join(
                    self.tmp_folder, "chunk-output-{0}.{1}".format(i, extension)
                )
            ]
        return cached_outputs


class BaseTabularFile(object):
    def __init__(
        self,
        path,
        IO,
        entity_is_list,
        expected_number,
        filter_by_column_validity=None,
        sniff_line_limit=100,
    ):
        self.logger = logger
        self.path = os.path.abspath(path)
        self.IO = IO
        self.entity_is_list = entity_is_list
        self.expected_number = expected_number
        self.filter_by_column_validity = filter_by_column_validity
        self.sniff_line_limit = sniff_line_limit
        self._has_header = None
        self._data = None
        self._string_delimiter = self.get_string_delimiter()
        self._column_delimiter = self.get_delimiter()
        self.logger.debug("Expected number: {0}".format(self.expected_number))
        self.logger.debug("Entity is list: {0}".format(self.entity_is_list))

    def _get_delimiter_by_extension(self):
        if self.path.endswith(".csv"):
            return ","
        if self.path.endswith(".tsv"):
            return "\t"
        return DEFAULT_DELIMITER

    def get_delimiter(self):
        delimiters = collections.defaultdict(int)
        default_extension = self._get_delimiter_by_extension()
        with open(self.path, "r") as f:
            i = 0
            R = []
            for l in f:
                R += [l]
                i += 1
                if i > self.sniff_line_limit:
                    self.logger.debug("Stopping sniffer for finding delimiter")
                    break
                try:
                    delimiter = csv.Sniffer().sniff(l, delimiters="\t,;").delimiter
                except:
                    delimiter = default_extension
                delimiters[delimiter] += 1
        delimiters = [
            (k, v) for k, v in sorted(delimiters.items(), key=lambda item: -item[1])
        ]
        return delimiters[0][0]

    def get_string_delimiter(self):
        return self.IO.string_delimiter()

    def is_key(self, v: str) -> bool:
        if not v:
            return False
        v = v.split(self._string_delimiter)
        n = len(v)
        t = 0
        for x in v:
            if self.IO.is_key(x):
                t += 1
        if t / n > MIN_COLUMN_VALIDITY:
            return True
        else:
            return False

    def is_input(self, v: str) -> bool:
        if not v:
            return False
        v = v.split(self._string_delimiter)
        n = len(v)
        t = 0
        for x in v:
            if self.IO.is_input(x):
                t += 1
        if t / n > MIN_COLUMN_VALIDITY:
            return True
        else:
            return False

    def resolve_columns(self):
        input = collections.defaultdict(int)
        key = collections.defaultdict(int)
        with open(self.path, "r") as f:
            reader = csv.reader(f, delimiter=self._column_delimiter)
            N = 0
            for i, r in enumerate(reader):
                if len(r) == 1:
                    self.matching = {"input": [0], "key": None}
                    return self.matching
                if i > self.sniff_line_limit:
                    self.logger.debug("Stopping sniffer for resolving column types")
                    break
                for j, v in enumerate(r):
                    if self.is_input(v):
                        input[j] += 1
                    if self.is_key(v):
                        key[j] += 1
                N += 1
        if len(key) == 0:
            key = None
        else:
            keys_and_counts = [
                (k, v)
                for k, v in sorted(key.items(), key=lambda item: -item[1])
                if v > 0
            ]
            keys_and_counts = keys_and_counts[: self.expected_number]
            key = []
            for k, n in keys_and_counts:
                if self.filter_by_column_validity:
                    if (float(n + 1) / N) < MIN_COLUMN_VALIDITY:
                        continue
                    else:
                        key += [k]
                else:
                    key += [k]
            if len(key) == 0:
                key = None
            else:
                key = sorted(key)
        if len(input) == 0:
            input = None
        else:
            inputs_and_counts = [
                (k, v)
                for k, v in sorted(input.items(), key=lambda item: -item[1])
                if v > 0
            ]
            inputs_and_counts = inputs_and_counts[: self.expected_number]
            input = []
            for k, n in inputs_and_counts:
                if self.filter_by_column_validity:
                    if (float(n + 1) / N) < MIN_COLUMN_VALIDITY:
                        continue
                    else:
                        input += [k]
                else:
                    input += [k]
            if len(input) == 0:
                input = None
            else:
                input = sorted(input)
        if key is not None:
            assert len(key) == self.expected_number
        if input is not None:
            assert len(input) == self.expected_number
        if key is not None and input is not None:
            assert len(key) == len(input)
        self.matching = {"input": input, "key": key}
        return self.matching

    def has_header(self):
        if self._has_header is not None:
            return self._has_header
        self.resolve_columns()
        with open(self.path, "r") as f:
            reader = csv.reader(f, delimiter=self._column_delimiter)
            candidate_header = next(reader)
        input = self.matching["input"]
        if input is not None:
            hi = []
            for i in input:
                v = candidate_header[i]
                hi += [self.is_input(v)]
            hi = np.any(hi)
        else:
            hi = False
        key = self.matching["key"]
        if key is not None:
            hk = []
            for k in key:
                v = candidate_header[k]
                hk += [self.is_key(v)]
            hk = np.any(hk)
        else:
            hk = False
        if hi or hk:
            self._has_header = False
        else:
            self._has_header = True
        return self._has_header

    def read_input_columns(self):
        if self._data is not None:
            return self._data
        header = self.has_header()
        self.logger.debug("Has header {0}".format(header))
        self.logger.debug("Schema {0}".format(self.matching))
        input = self.matching["input"]
        assert input is not None
        with open(self.path, "r") as f:
            R = []
            reader = csv.reader(f, delimiter=self._column_delimiter)
            if header:
                next(reader)
            for l in reader:
                r = []
                for i in input:
                    r += [l[i]]
                R += [r]
        self._data = R
        return self._data

    def is_single_input(self):
        if self._data is None:
            data = self.read_input_columns()
        else:
            data = self._data
        if self.entity_is_list:
            if len(data) == 1:
                return True
            else:
                if self.expected_number == 1:
                    for r in data:
                        r = r[0]
                        if self._string_delimiter in r:
                            return False
                    return True
                else:
                    if len(data[0]) > 1:
                        return True
                    else:
                        return False
        else:
            if len(data) == 1:
                return True
            else:
                return False

    def is_flattened(self):
        data = self.read_input_columns()
        flat_evidence = 0
        total_evidence = 0
        for r in data:
            for x in r:
                total_evidence += 1
                if self._string_delimiter in x:
                    flat_evidence += 1
        if (flat_evidence / total_evidence) > FLATTENED_EVIDENCE:
            return True
        else:
            return False


class TabularFileShapeStandardizer(BaseTabularFile):
    def __init__(self, src_path, dst_path, input_shape, IO, sniff_line_limit=100):
        if type(input_shape) is str:
            self.input_shape = InputShape(input_shape).get()
        else:
            self.input_shape = input_shape
        if type(self.input_shape) is InputShapeSingle:
            expected_number = 1
            entity_is_list = False
            filter_by_column_validity = True
            self._standardizer = self._standardize_single
        if type(self.input_shape) is InputShapeList:
            expected_number = 1
            entity_is_list = True
            filter_by_column_validity = True
            self._standardizer = self._standardize_list
        if type(self.input_shape) is InputShapePairOfLists:
            expected_number = 2
            entity_is_list = True
            filter_by_column_validity = False
            self._standardizer = self._standardize_pair_of_lists
        self.src_path = os.path.abspath(src_path)
        self.dst_path = os.path.abspath(dst_path)
        BaseTabularFile.__init__(
            self,
            path=self.src_path,
            IO=IO,
            entity_is_list=entity_is_list,
            expected_number=expected_number,
            filter_by_column_validity=filter_by_column_validity,
            sniff_line_limit=sniff_line_limit,
        )
        self.dst_string_delimiter = self._string_delimiter
        self.dst_column_delimiter = IO.column_delimiter()
        self.read_input_columns()

    def _standardize_single(self):
        self.logger.debug("Standardizing input single")
        with open(self.dst_path, "w") as f:
            writer = csv.writer(f, delimiter=self.dst_column_delimiter)
            writer.writerow(["input_0"])
            for r in self._data:
                writer.writerow(r)

    def _standardize_list(self):
        self.logger.debug("Standardizing input list")
        if self.is_single_input():
            self.logger.debug("This seems to be a single input!")
            R = []
            for r in self._data:
                R += r
            R = [[self.dst_string_delimiter.join(R)]]
        else:
            self.logger.debug("More than one input has been found")
            R = []
            for r in self._data:
                R += [r]
        with open(self.dst_path, "w") as f:
            writer = csv.writer(f, delimiter=self.dst_column_delimiter)
            writer.writerow(["input_0"])
            for r in R:
                writer.writerow(r)

    def _standardize_pair_of_lists(self):
        self.logger.debug("Standardizing input pair of lists")
        if self.is_single_input():
            self.logger.debug("This seems to be a single input!")
            if self.is_flattened():
                self.logger.debug("But it is flattened")
                R = []
                for r in self._data:
                    R += [r]
            else:
                self.logger.debug("Not flattened. Flattening now!")
                R_0 = []
                R_1 = []
                for r in self._data:
                    if r[0]:
                        R_0 += [r[0]]
                    if r[1]:
                        R_1 += [r[1]]
                R = [
                    [
                        self.dst_string_delimiter.join(R_0),
                        self.dst_string_delimiter.join(R_1),
                    ]
                ]
        else:
            self.logger.debug("More than one input has been found")

        with open(self.dst_path, "w") as f:
            self.logger.debug("Writing standardized input to {0}".format(self.dst_path))
            writer = csv.writer(f, delimiter=self.dst_column_delimiter)
            writer.writerow(["input_0", "input_1"])
            for r in R:
                writer.writerow(r)

    def standardize(self):
        if type(self.input_shape) is InputShapeSingle:
            self._standardize_single()
        if type(self.input_shape) is InputShapeList:
            self._standardize_list()
        if type(self.input_shape) is InputShapePairOfLists:
            self._standardize_pair_of_lists()


class StandardTabularFileReader(BatchCacher):
    def __init__(self, path):
        BatchCacher.__init__(self)
        self.path = os.path.abspath(path)
        self.logger = logger
        self.logger.debug("Reading standard file from {0}".format(self.path))
        self._column_delimiter = self.get_delimiter()
        self._has_header = True

    def get_delimiter(self):
        if self.path.endswith(".csv"):
            return ","
        if self.path.endswith(".tsv"):
            return "\t"
        return ","

    def read_header(self):
        with open(self.path, "r") as f:
            return csv.reader(f, delimiter=self._column_delimiter)

    def read(self):
        with open(self.path, "r") as f:
            R = []
            reader = csv.reader(f, delimiter=self._column_delimiter)
            next(reader)
            for r in reader:
                R += [r]
        return R

    def is_worth_splitting(self):
        with open(self.path, "r") as f:
            n = 0
            for _ in f:
                n += 1
        self.logger.debug("File has {0} lines".format(n))
        if n > FILE_CHUNKSIZE:
            self.logger.debug("Worth splitting it")
            return True
        else:
            return False

    def split_in_cache(self):
        ft = FileTyper(self.path)
        extension = ft.get_extension()
        self.logger.debug("Splitting file in cache: {0}".format(self.tmp_folder))
        has_header = self._has_header
        with open(self.path, "r") as f:
            j = 0
            if has_header:
                header = next(f)
            i = 0
            for l in f:
                if i == 0:
                    gn = os.path.join(
                        self.tmp_folder, "chunk-input-{0}.{1}".format(j, extension)
                    )
                    g = open(gn, "w")
                    if has_header:
                        g.write(header)
                    g.write(l)
                    i += 1
                elif i == FILE_CHUNKSIZE:
                    g.write(l)
                    g.close()
                    j += 1
                    i = 0
                else:
                    g.write(l)
                    i += 1
        return self.get_cached_input_files()


class TabularFileReader(StandardTabularFileReader):
    def __init__(self, path, IO, sniff_line_limit=100):
        self.src_path = os.path.abspath(path)
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        self.dst_path = os.path.join(self.tmp_folder, "standard_input_file.csv")
        self.path = self.dst_path
        self.IO = IO
        self.input_shape = IO.input_shape
        self.sniff_line_limit = sniff_line_limit
        self._string_delimiter = IO.string_delimiter()
        if type(self.input_shape) is InputShapeSingle:
            self._datum_parser = self._datum_parser_single
        if type(self.input_shape) is InputShapeList:
            self._datum_parser = self._datum_parser_list
        if type(self.input_shape) is InputShapePairOfLists:
            self._datum_parser = self._datum_parser_pair_of_lists
        self._standardize()
        StandardTabularFileReader.__init__(self, path=self.dst_path)

    def _standardize(self):
        tfss = TabularFileShapeStandardizer(
            src_path=self.src_path,
            dst_path=self.dst_path,
            input_shape=self.input_shape,
            IO=self.IO,
            sniff_line_limit=self.sniff_line_limit,
        )
        tfss.standardize()

    def _datum_parser_single(self, r):
        return r[0]

    def _datum_parser_list(self, r):
        return r[0].split(self._string_delimiter)

    def _datum_parser_pair_of_lists(self, r):
        return [r[0].split(self._string_delimiter), r[1].split(self._string_delimiter)]

    def read(self):
        if not os.path.exists(self.dst_path):
            self._standardize()
        with open(self.path, "r") as f:
            R = []
            reader = csv.reader(f, delimiter=self._column_delimiter)
            next(reader)
            for r in reader:
                R += [self._datum_parser(r)]
        return R


class BaseJsonFile(object):
    def __init__(self, path, IO, entity_is_list, expected_number):
        self.logger = logger
        self.path = os.path.abspath(path)
        self.logger.debug("Reading from JSON file {0}".format(self.path))
        self.entity_is_list = entity_is_list
        self.expected_number = expected_number
        self.IO = IO
        self.logger.debug("Expected number: {0}".format(self.expected_number))
        self.logger.debug("Entity is list: {0}".format(self.entity_is_list))
        self._data = None

    def read_input_json(self):
        if self._data is None:
            with open(self.path, "r") as f:
                self._data = json.load(f)
        return self._data

    def is_single_input(self):
        if self._data is None:
            data = self.read_input_json()
        else:
            data = self._data
        if self.entity_is_list:
            assert type(data) is list
            if self.expected_number == 1:
                one_element = data[0]
                if type(one_element) is str:
                    return True
                else:
                    return False
            else:
                one_element = data[0]
                assert type(one_element) is list
                one_inner_element = one_element[0]
                if type(one_inner_element) is list:
                    return False
                else:
                    return True
        else:
            if type(data) is str:
                return True
            else:
                return False


class JsonFileShapeStandardizer(BaseJsonFile):
    def __init__(self, src_path, dst_path, input_shape, IO):
        self.src_path = os.path.abspath(src_path)
        self.dst_path = os.path.abspath(dst_path)
        if type(input_shape) is str:
            self.input_shape = InputShape(input_shape).get()
        else:
            self.input_shape = input_shape
        if type(self.input_shape) is InputShapeSingle:
            expected_number = 1
            entity_is_list = False
        if type(self.input_shape) is InputShapeList:
            expected_number = 1
            entity_is_list = True
        if type(self.input_shape) is InputShapePairOfLists:
            expected_number = 2
            entity_is_list = True
        BaseJsonFile.__init__(
            self,
            path=self.src_path,
            IO=IO,
            entity_is_list=entity_is_list,
            expected_number=expected_number,
        )

    def standardize(self):
        if self.is_single_input():
            data = [self.read_input_json()]
        else:
            data = self.read_input_json()
        with open(self.dst_path, "w") as f:
            json.dump(data, f, indent=4)


class StandardJsonFileReader(BatchCacher):
    def __init__(self, path):
        BatchCacher.__init__(self)
        self.path = os.path.abspath(path)
        self.logger = logger
        self.logger.debug("Reading standard file from {0}".format(self.path))

    @staticmethod
    def _chunker(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    def read(self):
        with open(self.path, "r") as f:
            return json.load(f)

    def is_worth_splitting(self):
        n = len(self.read())
        self.logger.debug("File has {0} entris".format(n))
        if n > FILE_CHUNKSIZE:
            self.logger.debug("Worth splitting it")
            return True
        else:
            return False

    def split_in_cache(self):
        self.logger.debug("Splitting file in cache: {0}".format(self.tmp_folder))
        with open(self.path, "r") as f:
            data = json.load(f)
        for i, data_chunk in enumerate(self._chunker(data, FILE_CHUNKSIZE)):
            g = os.path.join(self.tmp_folder, "chunk-input-{0}.json".format(i))
            with open(g, "w") as f:
                json.dump(data_chunk, g, indent=4)
        return self.get_cached_input_files()


class JsonFileReader(StandardJsonFileReader):
    def __init__(self, path, IO):
        self.src_path = os.path.abspath(path)
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        self.dst_path = os.path.join(self.tmp_folder, "standard_input_file.json")
        self.path = self.dst_path
        self.IO = IO
        self.input_shape = IO.input_shape
        self._standardize()
        StandardJsonFileReader.__init__(self, path=self.dst_path)

    def _standardize(self):
        jfss = JsonFileShapeStandardizer(
            src_path=self.src_path,
            dst_path=self.dst_path,
            input_shape=self.input_shape,
            IO=self.IO,
        )
        jfss.standardize()

    def read(self):
        if not os.path.exists(self.dst_path):
            self._standardize()
        with open(self.path, "r") as f:
            return json.load(f)
