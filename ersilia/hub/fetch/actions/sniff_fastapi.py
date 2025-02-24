import csv
import json
import os
from pathlib import Path

from .... import ErsiliaBase, throw_ersilia_exception
from ....default import API_SCHEMA_FILE, MODEL_SIZE_FILE
from ....utils.exceptions_utils.fetch_exceptions import (
    SniffFastApiColumnsDontMatch,
    SniffFastApiColumnTypesIncompatibility,
)
from . import BaseAction

# TODO move this to default?
BUILTIN_EXAMPLE_INPUT_FILE_NAME = ["run_input.csv", "input.csv"]
BUILTIN_EXAMPLE_OUTPUT_FILE_NAME = ["run_output.csv", "output.csv"]
N = 3

# TODO move this to default
COLUMNS_FILENAME = ["run_columns.csv"]


class BuiltinExampleReader(ErsiliaBase):
    """
    Reads built-in examples for a FastAPI model.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    config_json : dict
        Configuration settings for the reader.

    Methods
    -------
    input() -> list
        Returns a list of input examples.
    output() -> list
        Returns a list of output examples.
    columns() -> list
        Returns a list of columns from the output file.
    """

    def __init__(self, model_id: str, config_json: dict):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.example_input_file = None
        self.example_output_file = None
        for example_input_filename in BUILTIN_EXAMPLE_INPUT_FILE_NAME:
            full_path = os.path.join(
                self._get_bundle_location(self.model_id),
                "model",
                "framework",
                "examples",
                example_input_filename,
            )
            if os.path.exists(full_path):
                self.example_input_file = full_path
        for example_output_filename in BUILTIN_EXAMPLE_OUTPUT_FILE_NAME:
            full_path = os.path.join(
                self._get_bundle_location(self.model_id),
                "model",
                "framework",
                "examples",
                example_output_filename,
            )
            if os.path.exists(full_path):
                self.example_output_file = full_path

    def input(self) -> list:
        """
        Returns a list of input examples.

        Returns
        -------
        list
            List of input examples.
        """
        if self.example_input_file is None:
            return None
        data = []
        with open(self.example_input_file, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for r in reader:
                data += [r[0]]
        return data[:N]

    def output(self) -> list:
        """
        Returns a list of output examples.

        Returns
        -------
        list
            List of output examples.
        """
        if self.example_output_file is None:
            return None
        data = []
        with open(self.example_output_file, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for r in reader:
                data += [r[0]]
        return data[:N]

    def columns(self) -> list:
        """
        Returns a list of columns.

        Returns
        -------
        list
            List of columns.
        """
        if self.example_output_file is None:
            return None
        with open(self.example_output_file, "r") as f:
            reader = csv.reader(f)
            columns = next(reader)
        return columns


class ColumnsFileReader(ErsiliaBase):
    """
    Reads column names

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    config_json : dict
        Configuration settings for the reader.

    Methods
    -------
    names() -> list
        Returns a list of column names.
    types() -> list
        Returns a list of column types.
    """

    def __init__(self, model_id: str, config_json: dict):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        full_path = os.path.join(
            self._get_bundle_location(self.model_id),
            "model",
            "framework",
            "columns",
            "run_columns.csv",
        )
        if os.path.exists(full_path):
            self.columns_file = os.path.join(full_path)
        else:
            self.columns_file = None

    def names(self) -> list:
        """
        Returns a list of column names.

        Returns
        -------
        list
            List of column names.
        """
        if self.columns_file is None:
            return None
        names = []
        with open(self.columns_file, "r") as f:
            reader = csv.reader(f)
            header = next(reader)
            idx = header.index("name")
            for r in reader:
                names += [r[idx]]
        return names

    def types(self) -> list:
        """
        Returns a list of column types.

        Returns
        -------
        list
            List of column types.
        """
        if self.columns_file is None:
            return None
        types = []
        with open(self.columns_file, "r") as f:
            reader = csv.reader(f)
            header = next(reader)
            idx = header.index("type")
            for r in reader:
                types += [r[idx]]
        return types


class ModelSniffer(BaseAction):
    """
    Infers the structure of a model by sniffing its inputs and outputs.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    config_json : dict
        Configuration settings for the sniffer.

    Methods
    -------
    sniff()
        Infers the structure of the model.
    """

    def __init__(self, model_id: str, config_json: dict):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )

    @staticmethod
    def _get_directory_size(dir: str) -> int:
        root_directory = Path(dir)
        bytes = sum(
            f.stat().st_size for f in root_directory.glob("**/*") if f.is_file()
        )
        return bytes

    def _get_size_in_mb(self) -> float:
        dest_dir = self._model_path(self.model_id)
        repo_dir = self._get_bundle_location(self.model_id)
        size = self._get_directory_size(dest_dir) + self._get_directory_size(repo_dir)
        mbytes = size / (1024**2)
        return mbytes

    @throw_ersilia_exception()
    def _get_schema_run(self) -> dict:
        er = BuiltinExampleReader(self.model_id, config_json=self.config_json)
        cr = ColumnsFileReader(self.model_id, config_json=self.config_json)
        names = cr.names()
        names_in_example = er.columns()
        if names != names_in_example:
            raise SniffFastApiColumnsDontMatch(model_id=self.model_id)
        types = set(cr.types())
        if types.issubset(set(["float", "integer"])):
            output_type = "numeric_array"
        elif types.issubset(set(["string"])):
            output_type = "string_array"
        else:
            raise SniffFastApiColumnTypesIncompatibility(model_id=self.model_id)
        data = {
            "input": {
                "key": {"type": "string"},
                "input": {"type": "string"},
                "text": {"type": "string"},
            },
            "output": {
                "outcome": {"type": output_type, "shape": [len(names)], "meta": names}
            },
        }
        return data

    def _get_schema(self) -> dict:
        schema = {"run": self._get_schema_run()}
        return schema

    def sniff(self):
        """
        Infers the structure of the model by sniffing by reading the columns file.

        This method:
        - Calculates and saves the model size.
        - Creates an API schema
        """
        self.logger.debug("Sniffing model")
        self.logger.debug("Getting model size")
        size = self._get_size_in_mb()
        self.logger.debug("Model size is {0} MB".format(size))
        path = os.path.join(self._model_path(self.model_id), MODEL_SIZE_FILE)
        with open(path, "w") as f:
            json.dump({"size": size, "units": "MB"}, f, indent=4)
        self.logger.debug("Resolving API schema from columns file...")
        all_schemas = self._get_schema()
        path = os.path.join(self._model_path(self.model_id), API_SCHEMA_FILE)
        with open(path, "w") as f:
            json.dump(all_schemas, f, indent=4)
        self.logger.debug("API schema saved at {0}".format(path))
