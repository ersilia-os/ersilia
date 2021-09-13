import os
import json
import importlib
import tempfile

from .readers.file import TabularFileReader
from .. import ErsiliaBase

ERSILIA_CFG = "ersilia.json"


class BaseIOGetter(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _read_ersilia_file(self, model_id):
        path = os.path.join(self._model_path(model_id), ERSILIA_CFG)
        if not os.path.exists(path):
            return None
        with open(path, "r") as f:
            return json.load(f)

    def _read_input_type(self, model_id):
        data = self._read_ersilia_file(model_id)
        if data is None:
            return None
        else:
            return data["input"]

    def get(self, model_id):
        input_type = self._read_input_type(model_id)
        if input_type is None:
            input_type = "naive"
        module = ".types.{0}".format(input_type)
        return importlib.import_module(module, package="ersilia.io").IO


class _GenericAdapter(object):
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


class GenericInputAdapter(object):
    def __init__(self, model_id, config_json=None):
        baseio = BaseIOGetter(config_json=config_json).get(model_id)
        self.adapter = _GenericAdapter(baseio)

    def adapt(self, inp):
        data = self.adapter.adapt(inp)
        return data
