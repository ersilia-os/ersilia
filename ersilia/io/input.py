import os
import json

from .readers.file import TabularFileReader


class GenericAdapter(object):
    def __init__(self, BaseIO):
        self.IO = BaseIO()

    def _is_file(self, inp):
        if os.path.isfile(inp):
            return True
        else:
            return False

    def _is_list(self, inp):
        if type(inp) is list:
            return True
        else:
            return False

    def _is_string(self, inp):
        if type(inp) is str:
            return True
        else:
            return False

    def _string_reader(self, inp):
        try:
            data = eval(inp)
        except:
            data = inp
        if type(data) is str:
            data = [data]
        return data

    def _file_reader(self, inp):
        reader = TabularFileReader(self.IO)
        data = reader.read(inp)
        return data

    def adapt(self, inp):
        """Given an input object (typically text), it returns a JSON Serializable object (typically a list of dictionaries)"""
        if self._is_string(inp):
            if self._is_file(inp):
                data = self._file_reader(inp)
            else:
                data = self._string_reader(inp)
        elif self._is_list(inp):
            data = inp
        else:
            return None
        data = [self.IO.parse(d) for d in data]
        return data
