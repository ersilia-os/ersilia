import os
import csv
import json
from pathlib import Path

from .... import throw_ersilia_exception

from . import BaseAction
from .... import ErsiliaBase
from ....default import MODEL_SIZE_FILE


N = 3

BUILTIN_EXAMPLE_FILE_NAME = "example.csv"
BUILTIN_OUTPUT_FILE_NAME = "output.csv"


# TODO for now this is not used
class BuiltinExampleReader(ErsiliaBase):
    def __init__(self, model_id, config_json):
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

    def input_example(self):
        data = []
        with open(self.input_example_file, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for r in reader:
                data += [r[0]]
        return data[:N]

    def output_example(self):
        data = []
        with open(self.input_example_file, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for r in reader:
                data += [r[0]]
        return data[:N]


class ModelSniffer(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )

    @staticmethod
    def _get_directory_size(dir):
        root_directory = Path(dir)
        bytes = sum(
            f.stat().st_size for f in root_directory.glob("**/*") if f.is_file()
        )
        return bytes

    def _get_size_in_mb(self):
        dest_dir = self._model_path(self.model_id)
        repo_dir = self._get_bundle_location(self.model_id)
        size = self._get_directory_size(dest_dir) + self._get_directory_size(repo_dir)
        mbytes = size / (1024**2)
        return mbytes

    @throw_ersilia_exception
    def sniff(self):
        self.logger.debug("Sniffing model")
        self.logger.debug("Getting model size")
        size = self._get_size_in_mb()
        self.logger.debug("Model size is {0} MB".format(size))
        path = os.path.join(self._model_path(self.model_id), MODEL_SIZE_FILE)
        with open(path, "w") as f:
            json.dump({"size": size, "units": "MB"}, f, indent=4)
