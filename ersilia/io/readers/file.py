# TODO: A file reader that automatically guesses a format based on the content of the file.

import csv


class TabularFileReader(object):
    def __init__(self, IO, sniff_line_limit=10):
        self.IO = IO
        self.sniff_line_limit = sniff_line_limit

    def get_dialect(self, path):
        with open(path, "r") as f:
            i = 0
            R = []
            for l in f:
                R += [l]
                i += 1
                if i > dialect_line_limit:
                    break
            self.dialect = csv.Sniffer().sniff(R, delimiters="\t,;")

    def has_header(self, path):
        pass

    def read(self, path):
        """Read a file.
        At the moment, this function only works for 1-column files having no header."""
        with open(path, "r") as f:
            R = []
            reader = csv.reader(f, delimiter="\t")
            for l in reader:
                R += [l[0]]
        return R
