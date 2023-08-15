import os
import json
import csv
import importlib
import itertools

from ..hub.content.card import ModelCard
from .. import ErsiliaBase
from .. import throw_ersilia_exception

from ..utils.exceptions_utils.exceptions import NullModelIdentifierError

from .shape import InputShape
from .shape import InputShapeSingle, InputShapeList, InputShapePairOfLists
from .readers.pyinput import PyInputReader
from .readers.file import TabularFileReader, JsonFileReader


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

    def _read_shape_from_card(self, model_id):
        self.logger.debug("Reading shape from {0}".format(model_id))
        try:
            input_shape = self.mc.get(model_id)["Input Shape"]
        except:
            input_shape = None
        self.logger.debug("Input Shape: {0}".format(input_shape))
        return InputShape(input_shape).get()

    def shape(self, model_id):
        return self._read_shape_from_card(model_id)

    def _get_from_model(self, model_id):
        input_type = self._read_input_from_card(model_id)
        if input_type is None:
            input_type = "naive"
        input_shape = self.shape(model_id)
        self.logger.debug("Input type is: {0}".format(input_type))
        self.logger.debug("Input shape is: {0}".format(input_shape.name))
        module = ".types.{0}".format(input_type)
        self.logger.debug("Importing module: {0}".format(module))
        return importlib.import_module(module, package="ersilia.io").IO(
            input_shape=input_shape
        )

    def _get_from_specs(self, input_type, input_shape):
        input_type = input_type.lower()
        input_shape = InputShape(input_shape).get()
        module = ".types.{0}".format(input_type)
        return importlib.import_module(module, package="ersilia.io").IO(
            input_shape=input_shape
        )

    def get(self, model_id=None, input_type=None, input_shape=None):
        if model_id is not None:
            return self._get_from_model(model_id=model_id)
        else:
            return self._get_from_specs(input_type=input_type, input_shape=input_shape)


class _GenericAdapter(object):
    def __init__(self, BaseIO):
        self.IO = BaseIO

    def _is_file(self, inp):
        if not self._is_string(inp):
            return False
        if os.path.isfile(inp):
            return True
        else:
            return False

    def _is_python_instance(self, inp):
        if type(inp) is str:
            if os.path.isfile(inp):
                return False
        return True

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

    def _try_to_eval(self, inp):
        try:
            data = eval(inp)
        except:
            data = inp
        return data

    def _is_tabular_file(self, inp):
        if inp.endswith(".csv") or inp.endswith(".tsv"):
            return True
        else:
            return False

    def _is_json_file(self, inp):
        if inp.endswith(".json"):
            return True
        else:
            return False

    def _py_input_reader(self, inp):
        reader = PyInputReader(input=inp, IO=self.IO)
        data = reader.read()
        return data

    def _file_reader(self, inp):
        reader = None
        if self._is_tabular_file(inp):
            reader = TabularFileReader(path=inp, IO=self.IO)
        if self._is_json_file(inp):
            reader = JsonFileReader(path=inp, IO=self.IO)
        data = reader.read()
        return data

    def _adapt(self, inp):
        if self._is_file(inp):
            return self._file_reader(inp)
        inp = self._try_to_eval(inp)
        if self._is_python_instance(inp):
            return self._py_input_reader(inp)

    def adapt(self, inp):
        data = self._adapt(inp)
        data = [self.IO.parse(d) for d in data]
        return data


class GenericInputAdapter(object):
    def __init__(
        self, model_id=None, input_type=None, input_shape=None, config_json=None
    ):
        baseio = BaseIOGetter(config_json=config_json).get(
            model_id=model_id, input_type=input_type, input_shape=input_shape
        )
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

    def adapt_one_by_one(self, inp):
        data = self.adapter.adapt(inp)
        for d in data:
            yield d


class ExampleGenerator(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        self.check_model_id(model_id)
        self.IO = BaseIOGetter(config_json=config_json).get(model_id)
        ErsiliaBase.__init__(self, config_json=config_json)
        self.input_shape = self.IO.input_shape
        self._string_delimiter = self.IO.string_delimiter()
        self._force_simple = True
        if type(self.input_shape) is InputShapeSingle:
            self._flatten = self._flatten_single
            self._force_simple = False
        if type(self.input_shape) is InputShapeList:
            self._flatten = self._flatten_list
        if type(self.input_shape) is InputShapePairOfLists:
            self._flatten = self._flatten_pair_of_lists

    @throw_ersilia_exception
    def check_model_id(self, model_id):
        if model_id is None:
            raise NullModelIdentifierError(model=model_id)

    @staticmethod
    def _get_delimiter(file_name):
        extension = file_name.split(".")[-1]
        if extension == "tsv":
            return "\t"
        else:
            return ","

    def _flatten_single(self, datum):
        return [datum]

    def _flatten_list(self, datum):
        return [self._string_delimiter.join(datum)]

    def _flatten_pair_of_lists(self, datum):
        return [
            self._string_delimiter.join(datum[0]),
            self._string_delimiter.join(datum[1]),
        ]

    def test(self):
        return self.IO.test()

    def example(self, n_samples, file_name, simple):
        if not self._force_simple:
            simple = simple
        else:
            simple = True
        if file_name is None:
            data = [v for v in self.IO.example(n_samples)]
            if simple:
                data = [d["input"] for d in data]
            return data
        else:
            extension = file_name.split(".")[-1]
            if extension == "json":
                with open(file_name, "w") as f:
                    data = [v for v in self.IO.example(n_samples)]
                    if simple:
                        data = [d["input"] for d in data]
                    json.dump(data, f, indent=4)
            else:
                delimiter = self._get_delimiter(file_name)
                with open(file_name, "w", newline="") as f:
                    writer = csv.writer(f, delimiter=delimiter)
                    if simple:
                        for v in self.IO.example(n_samples):
                            writer.writerow(self._flatten(v["input"]))
                    else:
                        writer.writerow(["key", "input", "text"])
                        for v in self.IO.example(n_samples):
                            writer.writerow([v["key"], v["input"], v["text"]])
