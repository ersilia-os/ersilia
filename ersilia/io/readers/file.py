import collections
import csv
import json
import os

from ... import logger
from ...default import HEADER_INDICATORS
from ...utils.logging import make_temp_dir
from ..shape import InputShape, InputShapeList, InputShapePairOfLists, InputShapeSingle

MIN_COLUMN_VALIDITY = 0.8
FLATTENED_EVIDENCE = 0.2
DEFAULT_DELIMITER = ","
FILE_CHUNKSIZE = 10000


class FileTyper(object):
    """
    Class to determine the type of a file based on its extension.

    Parameters
    ----------
    path : str
        Path to the file.
    """

    def __init__(self, path):
        self.path = os.path.join(path)

    def is_valid_input_file(self):
        """
        Check if the file is a valid input file.

        Returns
        -------
        bool
            True if the file is a valid input file, False otherwise.
        """
        if self.is_csv() or self.is_tsv() or self.is_json():
            return True
        else:
            return False

    def is_valid_output_file(self):
        """
        Check if the file is a valid output file.

        Returns
        -------
        bool
            True if the file is a valid output file, False otherwise.
        """
        if self.is_csv() or self.is_tsv() or self.is_json() or self.is_hdf5():
            return True
        else:
            return False

    def is_tabular(self):
        """
        Check if the file is a tabular file (CSV or TSV).

        Returns
        -------
        bool
            True if the file is a tabular file, False otherwise.
        """
        if self.is_csv() or self.is_tsv():
            return True
        else:
            return False

    def is_csv(self):
        """
        Check if the file is a CSV file.

        Returns
        -------
        bool
            True if the file is a CSV file, False otherwise.
        """
        if self.path.endswith(".csv"):
            return True
        else:
            return False

    def is_tsv(self):
        """
        Check if the file is a TSV file.

        Returns
        -------
        bool
            True if the file is a TSV file, False otherwise.
        """
        if self.path.endswith(".tsv"):
            return True
        else:
            return False

    def is_hdf5(self):
        """
        Check if the file is an HDF5 file.

        Returns
        -------
        bool
            True if the file is an HDF5 file, False otherwise.
        """
        if self.path.endswith(".h5"):
            return True
        else:
            return False

    def is_json(self):
        """
        Check if the file is a JSON file.

        Returns
        -------
        bool
            True if the file is a JSON file, False otherwise.
        """
        if self.path.endswith(".json"):
            return True
        else:
            return False

    def get_extension(self):
        """
        Get the file extension.

        Returns
        -------
        str
            The file extension.
        """
        return self.path.split(".")[-1]


class BatchCacher(object):
    """
    Class to handle caching of file batches.
    """

    def __init__(self):
        self.tmp_folder = make_temp_dir(prefix="ersilia-")

    def get_cached_files(self, prefix):
        """
        Get cached files with a specific prefix.

        Parameters
        ----------
        prefix : str
            The prefix of the cached files.

        Returns
        -------
        list
            List of cached files with the specified prefix.
        """
        idx2fn = {}
        for fn in os.listdir(self.tmp_folder):
            if fn.startswith(prefix):
                idx = int(fn.split("-")[-1].split(".")[0])
                idx2fn[idx] = os.path.join(self.tmp_folder, fn)
        idxs = sorted(idx2fn.keys())
        return [idx2fn[idx] for idx in idxs]

    def get_cached_input_files(self):
        """
        Get cached input files.

        Returns
        -------
        list
            List of cached input files.
        """
        return self.get_cached_files(prefix="chunk-input-")

    def get_cached_output_files(self):
        """
        Get cached output files.

        Returns
        -------
        list
            List of cached output files.
        """
        return self.get_cached_files(prefix="chunk-output-")

    def name_cached_output_files(self, cached_inputs, output_template):
        """
        Name cached output files based on cached input files and an output template.

        Parameters
        ----------
        cached_inputs : list
            List of cached input files.
        output_template : str
            Template for naming the output files.

        Returns
        -------
        list
            List of named cached output files.
        """
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
    """
    Base class for handling tabular files.

    Parameters
    ----------
    path : str
        Path to the file.
    IO : object
        IO handler object.
    entity_is_list : bool
        Whether the entity is a list.
    expected_number : int
        Expected number of columns.
    filter_by_column_validity : bool, optional
        Whether to filter by column validity.
    sniff_line_limit : int, optional
        Line limit for sniffing the file.
    """

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
        self.header_indicators = set(HEADER_INDICATORS)

    def _get_delimiter_by_extension(self):
        if self.path.endswith(".csv"):
            return ","
        if self.path.endswith(".tsv"):
            return "\t"
        return DEFAULT_DELIMITER

    def get_delimiter(self):
        """
        Get the delimiter used in the file.

        Returns
        -------
        str
            The delimiter used in the file.
        """
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
        """
        Get the string delimiter used in the file.

        Returns
        -------
        str
            The string delimiter.
        """
        return self.IO.string_delimiter()

    def is_key(self, v: str) -> bool:
        """
        Check if a value is a key.

        Parameters
        ----------
        v : str
            The value to check.

        Returns
        -------
        bool
            True if the value is a key, False otherwise.
        """
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
        """
        Check if a value is an input.

        Parameters
        ----------
        v : str
            The value to check.

        Returns
        -------
        bool
            True if the value is an input, False otherwise.
        """
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
        """
        Resolve the columns in the file to determine input and key columns.

        Returns
        -------
        None
        """
        input = collections.defaultdict(int)
        key = collections.defaultdict(int)
        with open(self.path, "r") as f:
            reader = csv.reader(f, delimiter=self._column_delimiter)
            N = 0
            for i, r in enumerate(reader):
                if len(r) == 1:
                    self.matching = {"input": [0], "key": None}
                    self.logger.debug(
                        "Number of columns seems to be 1: assuming input is the only column: {0}".format(
                            self.matching
                        )
                    )
                if len(r) > 2:
                    raise Exception(
                        "Too many columns in the input file. Maximum number of columns is 2 (input and key)."
                    )
                if i > self.sniff_line_limit:
                    self.logger.debug("Stopping sniffer for resolving column types")
                    break
                for j, v in enumerate(r):
                    if self.is_input(v):
                        input[j] += 1
                    if self.is_key(v):
                        key[j] += 1
                N += 1
            self.logger.debug("Done with sniffing the file")
            self.logger.debug("Input: {0}".format(dict(input)))
            self.logger.debug("Key: {0}".format(dict(key)))
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
            self.logger.debug("Key: {0}".format(key))
        if len(input) == 0:
            input = None
        else:
            inputs_and_counts = [
                (k, v)
                for k, v in sorted(input.items(), key=lambda item: -item[1])
                if v > 0
            ]
            if key is not None:
                inputs_and_counts = [
                    (k, n) for k, n in inputs_and_counts if k not in key
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
            self.logger.debug("Input: {0}".format(input))
        if key is not None:
            assert len(key) == self.expected_number
        if input is not None:
            assert len(input) == self.expected_number
        if key is not None and input is not None:
            assert len(key) == len(input)
        self.matching = {"input": input, "key": key}

    def has_header(self):
        """
        Check if the file has a header.

        Returns
        -------
        bool
            True if the file has a header, False otherwise.
        """
        if self._has_header is not None:
            self.logger.debug("Has header is not None")
            return self._has_header
        self.resolve_columns()
        with open(self.path, "r") as f:
            reader = csv.reader(f, delimiter=self._column_delimiter)
            candidate_header = next(reader)
            self.logger.debug("Candidate header is {0}".format(candidate_header))
        self._has_header = False
        for c in candidate_header:
            if c.lower() in self.header_indicators:
                self._has_header = True
        return self._has_header

    def read_input_columns(self):
        """
        Read the input columns from the file.

        Returns
        -------
        list
            List of input columns.
        """
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
                h = next(reader)
            else:
                h = None
            if h is not None and len(h) == 1:
                for l in reader:
                    l = self._column_delimiter.join(l)
                    R += [[l]]
            else:
                for l in reader:
                    r = []
                    for i in input:
                        r += [l[i]]
                    R += [r]
        self._data = R
        return self._data

    def is_single_input(self):
        """
        Check if the file has a single input.

        Returns
        -------
        bool
            True if the file has a single input, False otherwise.
        """
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
        """
        Check if the file is flattened.

        Returns
        -------
        bool
            True if the file is flattened, False otherwise.
        """
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
    """
    Class to standardize the shape of tabular files.

    Parameters
    ----------
    src_path : str
        Source path of the file.
    dst_path : str
        Destination path of the standardized file.
    input_shape : str or object
        Input shape specification.
    IO : object
        IO handler object.
    sniff_line_limit : int, optional
        Line limit for sniffing the file.

    Examples
    --------
    .. code-block:: python

        tfss = TabularFileShapeStandardizer(
            "data.csv",
            "standard_data.csv",
            "single",
            IOHandler(),
        )
        tfss.standardize()
    """

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
            self.logger.debug("Writing standardized input to {0}".format(self.dst_path))
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
        """
        Standardize the shape of the tabular file.

        Returns
        -------
        None
        """
        if type(self.input_shape) is InputShapeSingle:
            self._standardize_single()
        if type(self.input_shape) is InputShapeList:
            self._standardize_list()
        if type(self.input_shape) is InputShapePairOfLists:
            self._standardize_pair_of_lists()


class StandardTabularFileReader(BatchCacher):
    """
    Class to read standard tabular files.

    Parameters
    ----------
    path : str
        Path to the file.
    """

    def __init__(self, path):
        BatchCacher.__init__(self)
        self.path = os.path.abspath(path)
        self.logger = logger
        self.logger.debug("Reading standard file from {0}".format(self.path))
        self._column_delimiter = self.get_delimiter()
        self._has_header = True

    def get_delimiter(self):
        """
        Get the delimiter used in the file.

        Returns
        -------
        str
            The delimiter used in the file.
        """
        if self.path.endswith(".csv"):
            return ","
        if self.path.endswith(".tsv"):
            return "\t"
        return ","

    def read_header(self):
        """
        Read the header of the file.

        Returns
        -------
        list
            List of header columns.
        """
        with open(self.path, "r") as f:
            return csv.reader(f, delimiter=self._column_delimiter)

    def read(self):
        """
        Read the content of the file.

        Returns
        -------
        list
            List of rows in the file.
        """
        with open(self.path, "r") as f:
            R = []
            reader = csv.reader(f, delimiter=self._column_delimiter)
            next(reader)
            for r in reader:
                R += [r]
        return R

    def is_worth_splitting(self):
        """
        Check if the file is worth splitting into smaller chunks.

        Returns
        -------
        bool
            True if the file is worth splitting, False otherwise.
        """
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
        """
        Split the file into smaller chunks and cache them.

        Returns
        -------
        list
            List of cached input files.
        """
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
    """
    Class to read and standardize tabular files.

    Parameters
    ----------
    path : str
        Path to the file.
    IO : object
        IO handler object.
    sniff_line_limit : int, optional
        Line limit for sniffing the file.
    """

    def __init__(self, path, IO, sniff_line_limit=100):
        self.src_path = os.path.abspath(path)
        self.tmp_folder = make_temp_dir(prefix="ersilia-")
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
        """
        Read the content of the file.

        Returns
        -------
        list
            List of rows in the file.
        """
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
    """
    Base class for handling JSON files.

    Parameters
    ----------
    path : str
        Path to the file.
    IO : object
        IO handler object.
    entity_is_list : bool
        Whether the entity is a list.
    expected_number : int
        Expected number of elements.
    """

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
        """
        Read the input JSON file.

        Returns
        -------
        dict or list
            Parsed JSON data.
        """
        if self._data is None:
            with open(self.path, "r") as f:
                self._data = json.load(f)
        return self._data

    def is_single_input(self):
        """
        Check if the JSON file has a single input.

        Returns
        -------
        bool
            True if the JSON file has a single input, False otherwise.
        """
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
    """
    Class to standardize the shape of JSON files.

    Parameters
    ----------
    src_path : str
        Source path of the file.
    dst_path : str
        Destination path of the standardized file.
    input_shape : str or object
        Input shape specification.
    IO : object
        IO handler object.
    """

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
        """
        Standardize the shape of the JSON file.

        Returns
        -------
        None
        """
        if self.is_single_input():
            data = [self.read_input_json()]
        else:
            data = self.read_input_json()
        with open(self.dst_path, "w") as f:
            json.dump(data, f, indent=4)


class StandardJsonFileReader(BatchCacher):
    """
    Class to read standard JSON files.

    Parameters
    ----------
    path : str
        Path to the file.

    Examples
    --------
    >>> sjfr = StandardJsonFileReader("data.json")
    >>> sjfr.read()
    [{'key': 'value'}, {'key': 'value'}]
    """

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
        """
        Read the content of the JSON file.

        Returns
        -------
        dict or list
            Parsed JSON data.
        """
        with open(self.path, "r") as f:
            return json.load(f)

    def is_worth_splitting(self):
        """
        Check if the JSON file is worth splitting into smaller chunks.

        Returns
        -------
        bool
            True if the JSON file is worth splitting, False otherwise.
        """
        n = len(self.read())
        self.logger.debug("File has {0} entris".format(n))
        if n > FILE_CHUNKSIZE:
            self.logger.debug("Worth splitting it")
            return True
        else:
            return False

    def split_in_cache(self):
        """
        Split the JSON file into smaller chunks and cache them.

        Returns
        -------
        list
            List of cached input files.
        """
        self.logger.debug("Splitting file in cache: {0}".format(self.tmp_folder))
        with open(self.path, "r") as f:
            data = json.load(f)
        for i, data_chunk in enumerate(self._chunker(data, FILE_CHUNKSIZE)):
            g = os.path.join(self.tmp_folder, "chunk-input-{0}.json".format(i))
            with open(g, "w") as f:
                json.dump(data_chunk, g, indent=4)
        return self.get_cached_input_files()


class JsonFileReader(StandardJsonFileReader):
    """
    Class to read and standardize JSON files.

    Parameters
    ----------
    path : str
        Path to the file.
    IO : object
        IO handler object.
    """

    def __init__(self, path, IO):
        self.src_path = os.path.abspath(path)
        self.tmp_folder = make_temp_dir(prefix="ersilia-")
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
        """
        Read the content of the JSON file.

        Returns
        -------
        dict or list
            Parsed JSON data.
        """
        if not os.path.exists(self.dst_path):
            self._standardize()
        with open(self.path, "r") as f:
            return json.load(f)
