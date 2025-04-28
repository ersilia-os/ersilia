import csv
import importlib
import itertools
import json
import os
import shutil

from click import secho

from .. import ErsiliaBase, throw_ersilia_exception
from ..default import PREDEFINED_EXAMPLE_FILES
from ..hub.content.card import ModelCard
from ..utils.exceptions_utils.exceptions import NullModelIdentifierError
from .readers.file import JsonFileReader, TabularFileReader
from .readers.pyinput import PyInputReader
from .shape import InputShape, InputShapeList, InputShapePairOfLists, InputShapeSingle


class BaseIOGetter(ErsiliaBase):
    """
    Base class to get IO handlers based on model or specifications.

    Parameters
    ----------
    config_json : dict, optional
        Configuration JSON.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.mc = ModelCard(config_json=config_json)

    def shape(self, model_id):
        """
        Get the input shape from the model card.

        Parameters
        ----------
        model_id : str
            Model identifier.

        Returns
        -------
        InputShape
            Input shape object.
        """
        return self._read_shape_from_card(model_id)

    def get(self, model_id=None, input_type=None, input_shape=None):
        """
        Get the IO handler based on model or specifications.

        Parameters
        ----------
        model_id : str, optional
            Model identifier.
        input_type : str, optional
            Input type.
        input_shape : str, optional
            Input shape.

        Returns
        -------
        IO
            IO handler object.
        """
        if model_id is not None:
            return self._get_from_model(model_id=model_id)
        else:
            return self._get_from_specs(input_type=input_type, input_shape=input_shape)

    def _read_input_from_card(self, model_id):
        self.logger.debug("Reading card from {0}".format(model_id))
        # This is because ersilia-pack adds another level in the JSON with the key "card"
        if "Input" not in self.mc.get(model_id):
            input_type = self.mc.get(model_id)["card"]["Input"]
        else:
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


class _GenericAdapter(object):
    """
    Class to adapt various input formats to a standard format.

    This class handles different types of inputs such as files, lists, and strings,
    and converts them into a standard format that can be processed by the IO handler.

    Parameters
    ----------
    BaseIO : object
        Base IO handler object.
    """

    def __init__(self, BaseIO):
        self.IO = BaseIO

    def adapt(self, inp):
        """
        Adapt the input to a standard format.

        Parameters
        ----------
        inp : any
            The input data.

        Returns
        -------
        list
            List of adapted data.
        """
        data = self._adapt(inp)
        data = [self.IO.parse(d) for d in data]
        return data

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


class GenericInputAdapter(object):
    """
    Class to adapt generic inputs to a standard format.

    This class uses the _GenericAdapter to handle different types of inputs and
    convert them into a standard format that can be processed in batches or one by one.

    Parameters
    ----------
    model_id : str, optional
        Model identifier.
    input_type : str, optional
        Input type.
    input_shape : str, optional
        Input shape.
    config_json : dict, optional
        Configuration JSON.
    """

    def __init__(
        self, model_id=None, input_type=None, input_shape=None, config_json=None
    ):
        baseio = BaseIOGetter(config_json=config_json).get(
            model_id=model_id, input_type=input_type, input_shape=input_shape
        )
        self.adapter = _GenericAdapter(baseio)

    def adapt(self, inp, batch_size):
        """
        Adapt the input data in batches.

        Parameters
        ----------
        inp : any
            The input data.
        batch_size : int
            Size of each batch.

        Yields
        ------
        list
            List of adapted data in batches.
        """
        data = self.adapter.adapt(inp)
        for chunk in self.batch_iter(data, batch_size):
            yield chunk

    def adapt_one_by_one(self, inp):
        """
        Adapt the input data one by one.

        Parameters
        ----------
        inp : any
            The input data.

        Yields
        ------
        dict
            Adapted data.
        """
        data = self.adapter.adapt(inp)
        for d in data:
            yield d

    def batch_iter(self, data, batch_size):
        """
        Yield batches of data.

        Parameters
        ----------
        data : iterable
            The data to be batched.
        batch_size : int
            The size of each batch.

        Yields
        ------
        iterable
            A batch of data.
        """
        it = iter(data)
        while True:
            chunk = tuple(itertools.islice(it, batch_size))
            if not chunk:
                break
            yield chunk


class ExampleGenerator(ErsiliaBase):
    """
    Class to generate examples for a model.

    Parameters
    ----------
    model_id : str
        Model identifier.
    config_json : dict, optional
        Configuration JSON.
    """

    def __init__(self, model_id, config_json=None):
        self.model_id = model_id
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

    def test(self):
        """
        Get test examples.

        Returns
        -------
        list
            List of test examples.
        """
        return self.IO.test()

    def random_example(self, n_samples, file_name, simple):
        """
        Generate random example data.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate.
        file_name : str
            File name to save the examples.
        simple : bool
            Whether to generate simple examples.

        Returns
        -------
        list or None
            List of example data or None if saved to file.
        """
        if not self._force_simple:
            simple = simple
        else:
            simple = True
        if file_name is None:
            data = [v for v in self.IO.example(n_samples)]
            if simple:
                data = [{"input": d["input"]} for d in data]
            return data
        else:
            extension = file_name.split(".")[-1]
            if extension == "json":
                with open(file_name, "w") as f:
                    data = [v for v in self.IO.example(n_samples)]
                    if simple:
                        data = [{"input": d["input"]} for d in data]
                    json.dump(data, f, indent=4)
            else:
                delimiter = self._get_delimiter(file_name)
                with open(file_name, "w", newline="") as f:
                    writer = csv.writer(f, delimiter=delimiter)
                    if simple:
                        writer.writerow(["input"])
                        for v in self.IO.example(n_samples):
                            writer.writerow(self._flatten(v["input"]))
                    else:
                        writer.writerow(["key", "input", "text"])
                        for v in self.IO.example(n_samples):
                            writer.writerow([v["key"], v["input"], v["text"]])

    def fixed_example(self, n_samples, file_name, simple):
        """
        Generate deterministic example data. Samples same types of examples for each any sampling.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate.
        file_name : str
            File name to save the examples.
        simple : bool
            Whether to generate simple examples.

        Returns
        -------
        list or None
            List of example data or None if saved to file.
        """
        if not self._force_simple:
            simple = simple
        else:
            simple = True
        if file_name is None:
            data = [v for v in self.IO.example_fixed(n_samples)]
            if simple:
                data = [{"input": d["input"]} for d in data]
            return data
        else:
            extension = file_name.split(".")[-1]
            if extension == "json":
                with open(file_name, "w") as f:
                    data = [v for v in self.IO.example_fixed(n_samples)]
                    if simple:
                        data = [{"input": d["input"]} for d in data]
                    json.dump(data, f, indent=4)
            else:
                delimiter = self._get_delimiter(file_name)
                with open(file_name, "w", newline="") as f:
                    writer = csv.writer(f, delimiter=delimiter)
                    if simple:
                        writer.writerow(["input"])
                        for v in self.IO.example_fixed(n_samples):
                            writer.writerow(self._flatten(v["input"]))
                    else:
                        writer.writerow(["key", "input", "text"])
                        for v in self.IO.example_fixed(n_samples):
                            writer.writerow([v["key"], v["input"], v["text"]])

    def predefined_example(self, file_name):
        """
        Get predefined example data.

        Parameters
        ----------
        file_name : str
            File name to save the examples.

        Returns
        -------
        bool
            True if predefined examples are available, False otherwise.
        """
        dest_folder = self._model_path(self.model_id)
        for pf in PREDEFINED_EXAMPLE_FILES:
            example_file = os.path.join(dest_folder, pf)
            if os.path.exists(example_file):
                shutil.copy(example_file, file_name)
                return True
            else:
                return False

    def example(self, n_samples, file_name, simple, try_predefined, deterministic):
        """
        Generate example data.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate.
        file_name : str
            File name to save the examples.
        simple : bool
            Whether to generate simple examples.
        try_predefined : bool
            Whether to try predefined examples first.

        Returns
        -------
        list or str
            List of example data or file content if saved to file.
        """
        predefined_available = False
        if try_predefined is True and file_name is not None:
            self.logger.debug("Trying with predefined input")
            predefined_available = self.predefined_example(file_name)

        if predefined_available:
            with open(file_name, "r") as f:
                return f.read()
        elif deterministic:
            secho(
                "Sampling input not randomly but in deterministic manner.",
                fg="green",
            )
            self.logger.debug("Sampling input not randomly but in deterministic manner")
            return self.fixed_example(
                n_samples=n_samples, file_name=file_name, simple=simple
            )
        else:
            if try_predefined and not predefined_available:
                secho(
                    "No predefined examples found for the model. Generating random examples.",
                    fg="yellow",
                )
            self.logger.debug("Randomly sampling input")
            return self.random_example(
                n_samples=n_samples, file_name=file_name, simple=simple
            )

    @throw_ersilia_exception()
    def check_model_id(self, model_id):
        """
        Check if the model ID is valid.

        Parameters
        ----------
        model_id : str
            Model identifier.

        Returns
        -------
        None

        Raises
        ------
        NullModelIdentifierError
            If the model ID is None.
        """
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
