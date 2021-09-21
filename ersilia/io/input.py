import os
import json
import importlib
import tempfile
import itertools

from .readers.file import TabularFileReader
from ..hub.content.card import ModelCard
from .. import ErsiliaBase


class BaseIOGetter(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.mc = ModelCard(config_json=config_json)

    def _read_input_from_card(self, model_id):
        self.logger.debug("Reading card from {0}".format(model_id))
        input_type = self.mc.get(model_id)["Input"]
        if len(input_type) != 1:
            self.logger.error("Ersilia does not deal with multiple inputs yet..!")
        else:
            input_type = input_type[0]
            return input_type.lower()

    def get(self, model_id):
        input_type = self._read_input_from_card(model_id)
        if input_type is None:
            input_type = "naive"
        self.logger.debug("Input type is: {0}".format(input_type))
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

    def batch_iter(self, data, batch_size):
        it = iter(data)
        while True:
            chunk = tuple(itertools.islice(it, batch_size))
            if not chunk:
                break
            yield chunk

    def adapt(self, inp, batch_size):
        data = self.adapter.adapt(inp)
        for chunk in self.batch_iter(data, batch_size):
            yield chunk
