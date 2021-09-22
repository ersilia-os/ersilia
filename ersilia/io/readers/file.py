import csv
import collections
from ... import logger

MIN_COLUMN_VALIDITY = 0.8


class TabularFileReader(object):
    def __init__(self, IO, sniff_line_limit=100):
        self.logger = logger
        self.IO = IO
        self.sniff_line_limit = sniff_line_limit

    def get_delimiter(self, path):
        delimiters = collections.defaultdict(int)
        with open(path, "r") as f:
            i = 0
            R = []
            for l in f:
                R += [l]
                i += 1
                if i > self.sniff_line_limit:
                    break
                delimiter = csv.Sniffer().sniff(l, delimiters="\t,;").delimiter
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
            for r in reader:
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
            return False
        else:
            return True

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
