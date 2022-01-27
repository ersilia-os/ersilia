import os
import tempfile
import csv
import collections
from ... import logger

MIN_COLUMN_VALIDITY = 0.8
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


class TabularFileReader(object):
    def __init__(self, IO, sniff_line_limit=100):
        self.logger = logger
        self.IO = IO
        self.sniff_line_limit = sniff_line_limit
        self._has_header = None

    def get_delimiter(self, path):
        delimiters = collections.defaultdict(int)
        with open(path, "r") as f:
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
                    delimiter = DEFAULT_DELIMITER
                delimiters[delimiter] += 1
        delimiters = [
            (k, v) for k, v in sorted(delimiters.items(), key=lambda item: -item[1])
        ]
        self.delimiter = delimiters[0][0]

    def resolve_columns(self, path):
        self.get_delimiter(path)
        input = collections.defaultdict(int)
        key = collections.defaultdict(int)
        with open(path, "r") as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            N = 0
            for i, r in enumerate(reader):
                if i > self.sniff_line_limit:
                    self.logger.debug("Stopping sniffer for resolving column types")
                    break
                for j, v in enumerate(r):
                    if self.IO.is_input(v):
                        input[j] += 1
                    if self.IO.is_key(v):
                        key[j] += 1
                N += 1
        if len(key) == 0:
            key = None
        else:
            key, n = [
                (k, v) for k, v in sorted(key.items(), key=lambda item: -item[1])
            ][0]
            if float(n + 1) / N < MIN_COLUMN_VALIDITY:
                key = None
        if len(input) == 0:
            input = None
        else:
            input, n = [
                (k, v) for k, v in sorted(input.items(), key=lambda item: -item[1])
            ][0]
            if float(n + 1) / N < MIN_COLUMN_VALIDITY:
                input = None
        matching = {"input": input, "key": key}
        self.matching = matching

    def has_header(self, path):
        if self._has_header is not None:
            return self._has_header
        self.resolve_columns(path)
        with open(path, "r") as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            candidate_header = next(reader)
        input = self.matching["input"]
        if input is not None:
            v = candidate_header[input]
            hi = self.IO.is_input(v)
        else:
            hi = False
        key = self.matching["key"]
        if key is not None:
            v = candidate_header[key]
            hk = self.IO.is_key(v)
        else:
            hk = False
        if hi or hk:
            self._has_header = False
        else:
            self._has_header = True
        return self._has_header

    def read(self, path):
        header = self.has_header(path)
        self.logger.debug("Has header {0}".format(header))
        self.logger.debug("Schema {0}".format(self.matching))
        input = self.matching["input"]
        with open(path, "r") as f:
            R = []
            reader = csv.reader(f, delimiter=self.delimiter)
            if header:
                next(reader)
            for l in reader:
                R += [l[input]]
        return R

    def is_worth_splitting(self, path):
        with open(path, "r") as f:
            n = 0
            for _ in f:
                n += 1
        self.logger.debug("File has {0} lines".format(n))
        if n > FILE_CHUNKSIZE:
            self.logger.debug("Worth splitting it")
            return True
        else:
            return False

    def split_in_cache(self, path):
        ft = FileTyper(path)
        extension = ft.get_extension()
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        self.logger.debug("Splitting file in cache: {0}".format(self.tmp_folder))
        has_header = self.has_header(path)
        with open(path, "r") as f:
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
