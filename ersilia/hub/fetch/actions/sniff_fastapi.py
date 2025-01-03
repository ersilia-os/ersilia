import csv
import json
import os
from pathlib import Path

from .... import ErsiliaBase, throw_ersilia_exception
from ....default import MODEL_SIZE_FILE
from . import BaseAction

N = 3

BUILTIN_EXAMPLE_FILE_NAME = "example.csv"
BUILTIN_OUTPUT_FILE_NAME = "output.csv"


# TODO for now this is not used
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
    input_example() -> list
        Returns a list of input examples.
    output_example() -> list
        Returns a list of output examples.
    """

    def __init__(self, model_id: str, config_json: dict):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.input_example_file = os.path.join(
            self._get_bundle_location(self.model_id),
            "model",
            "framework",
            BUILTIN_EXAMPLE_FILE_NAME,
        )
        self.output_example_file = os.path.join(
            self._get_bundle_location(self.model_id),
            "model",
            "framework",
            BUILTIN_OUTPUT_FILE_NAME,
        )

    def input_example(self) -> list:
        """
        Returns a list of input examples.

        Returns
        -------
        list
            List of input examples.
        """
        data = []
        with open(self.input_example_file, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for r in reader:
                data += [r[0]]
        return data[:N]

    def output_example(self) -> list:
        """
        Returns a list of output examples.

        Returns
        -------
        list
            List of output examples.
        """
        data = []
        with open(self.input_example_file, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for r in reader:
                data += [r[0]]
        return data[:N]


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
    def sniff(self):
        """
        Infers the structure of the model by sniffing its inputs and outputs.
        """
        self.logger.debug("Sniffing model")
        self.logger.debug("Getting model size")
        size = self._get_size_in_mb()
        self.logger.debug("Model size is {0} MB".format(size))
        path = os.path.join(self._model_path(self.model_id), MODEL_SIZE_FILE)
        with open(path, "w") as f:
            json.dump({"size": size, "units": "MB"}, f, indent=4)
