import os
import json
import tempfile
import types
import collections
import importlib
import __main__ as main
from .. import logger
from .base import ErsiliaBase
from .modelbase import ModelBase
from ..serve.autoservice import AutoService
from ..serve.schema import ApiSchema
from ..serve.api import Api
from ..io.input import ExampleGenerator
from ..default import MODEL_SIZE_FILE, CARD_FILE
from ..default import DEFAULT_BATCH_SIZE
from ..utils import tmp_pid_file
from ..utils.hdf5 import Hdf5DataLoader
from ..utils.terminal import yes_no_input
from ..lake.base import LakeBase

try:
    import pandas as pd
except:
    pd = None


class ErsiliaModel(ErsiliaBase):
    def __init__(
        self,
        model,
        save_to_lake=True,
        config_json=None,
        credentials_json=None,
        verbose=None,
        fetch_if_not_available=True,
    ):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.logger = logger
        if verbose is not None:
            if verbose:
                self.logger.set_verbosity(1)
            else:
                self.logger.set_verbosity(0)
        else:
            if not hasattr(main, "__file__"):
                self.logger.set_verbosity(0)
        self.save_to_lake = save_to_lake
        if self.save_to_lake:
            lake = LakeBase(config_json=self.config_json)
            if not lake.is_installed():
                self.logger.error(
                    "Isaura is not installed! Calculations will be done without storing and reading from the lake, unfortunately."
                )
                self.save_to_lake = False
        mdl = ModelBase(model)
        self._is_valid = mdl.is_valid()
        assert (
            self._is_valid
        ), "The identifier {0} is not valid. Please visit the Ersilia Model Hub for valid identifiers".format(
            model
        )
        self.config_json = config_json
        self.model_id = mdl.model_id
        self.slug = mdl.slug
        self.text = mdl.text
        self._is_available_locally = mdl.is_available_locally()
        if not self._is_available_locally and fetch_if_not_available:
            self.logger.info("Model is not available locally")
            do_fetch = yes_no_input(
                "Requested model {0} if not available locally. Do you want to fetch it? [Y/n]".format(
                    self.model_id
                ),
                default_answer="Y",
            )
            if do_fetch:
                fetch = importlib.import_module("ersilia.hub.fetch.fetch")
                mf = fetch.ModelFetcher(
                    config_json=self.config_json, credentials_json=self.credentials_json
                )
                mf.fetch(self.model_id)
            else:
                return
        self.api_schema = ApiSchema(
            model_id=self.model_id, config_json=self.config_json
        )
        self.autoservice = AutoService(
            model_id=self.model_id, config_json=self.config_json
        )
        self._set_apis()

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def is_valid(self):
        return self._is_valid

    def _set_api(self, api_name):
        def _method(input=None, output=None, batch_size=DEFAULT_BATCH_SIZE):
            return self.api(api_name, input, output, batch_size)

        setattr(self, api_name, _method)

    def _set_apis(self):
        apis_list = os.path.join(
            self._get_bundle_location(self.model_id), "apis_list.txt"
        )
        if os.path.exists(apis_list):
            with open(apis_list, "r") as f:
                for l in f:
                    api_name = l.rstrip()
                    self._set_api(api_name)
        else:
            with open(apis_list, "w") as f:
                for api_name in self.autoservice.service._get_apis_from_bento():
                    self._set_api(api_name)
                    f.write(api_name + os.linesep)
        self.apis_list = apis_list

    def _get_api_instance(self, api_name):
        model_id = self.model_id
        tmp_file = tmp_pid_file(model_id)
        assert os.path.exists(
            tmp_file
        ), "Process ID file does not exist. Please serve the model first!"
        with open(tmp_file, "r") as f:
            for l in f:
                url = l.rstrip().split()[1]
        if api_name is None:
            api_names = self.autoservice.get_apis()
            assert (
                len(api_names) == 1
            ), "More than one API found, please specificy api_name"
            api_name = api_names[0]
        api = Api(
            model_id=model_id,
            url=url,
            api_name=api_name,
            save_to_lake=self.save_to_lake,
            config_json=self.config_json,
        )
        return api

    def _api_runner_iter(self, api, input, output, batch_size):
        for result in api.post(input=input, output=output, batch_size=batch_size):
            assert (
                result is not None
            ), "Something went wrong. Please contact us at hello@ersila.io"
            yield result

    def _api_runner_return(self, api, input, output, batch_size):
        if output == "pandas" and pd is None:
            raise Exception
        if output == "json":
            R = []
            for r in self._api_runner_iter(
                api=api, input=input, output=None, batch_size=batch_size
            ):
                R += [r]
            return json.dumps(R, indent=4)
        else:
            tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
            tmp_output = os.path.join(tmp_folder, "temporary.h5")
            for r in self._api_runner_iter(
                api=api, input=input, output=tmp_output, batch_size=batch_size
            ):
                continue
            data = Hdf5DataLoader()
            data.load(tmp_output)
            if output == "numpy":
                return data.values[:]
            if output == "pandas":
                d = collections.OrderedDict()
                d["key"] = data.keys
                d["input"] = data.inputs
                for j, f in enumerate(data.features):
                    d[f] = data.values[:, j]
                return pd.DataFrame(d)
            if output == "dict":
                d = collections.OrderedDict()
                d["keys"] = data.keys
                d["inputs"] = data.inputs
                d["values"] = data.values
                return d

    @staticmethod
    def __output_is_file(output):
        if output is None:
            return False
        if type(output) != str:
            return False
        if output[-4:] == ".csv":
            return True
        if output[-3:] == ".h5":
            return True
        return False

    @staticmethod
    def __output_is_format(output):
        if output is None:
            return False
        if type(output) != str:
            return False
        if output == "json":
            return True
        if output == "numpy":
            return True
        if output == "pandas":
            return True
        if output == "dict":
            return True
        return False

    def _get_api_runner(self, output):
        if output is None:
            use_iter = True
        elif self.__output_is_file(output):
            use_iter = True
        elif self.__output_is_format(output):
            use_iter = False
        if use_iter:
            return self._api_runner_iter
        else:
            return self._api_runner_return

    def api(
        self, api_name=None, input=None, output=None, batch_size=DEFAULT_BATCH_SIZE
    ):
        api_instance = self._get_api_instance(api_name=api_name)
        api_runner = self._get_api_runner(output=output)
        result = api_runner(
            api=api_instance, input=input, output=output, batch_size=batch_size
        )
        if output is None:
            # Result is a generator
            return result
        if isinstance(result, types.GeneratorType):
            # Result is a file
            for r in result:
                pass
            return output
        else:
            # Result is a dict, a numpy array, a dataframe...
            return result

    def serve(self):
        self.autoservice.serve()
        self.url = self.autoservice.service.url
        self.pid = self.autoservice.service.pid

    def close(self):
        self.autoservice.close()

    def get_apis(self):
        return self.autoservice.get_apis()

    @property
    def paths(self):
        p = {
            "dest": self._model_path(self.model_id),
            "repository": self._get_bundle_location(self.model_id),
            "bentoml": self._get_bentoml_location(self.model_id),
        }
        return p

    @property
    def input_type(self):
        with open(os.path.join(self._model_path(self.model_id), CARD_FILE), "r") as f:
            return [x.lower() for x in json.load(f)["Input"]]

    @property
    def output_type(self):
        with open(os.path.join(self._model_path(self.model_id), CARD_FILE), "r") as f:
            return [x.lower() for x in json.load(f)["Output"]]

    @property
    def schema(self):
        return self.api_schema.schema

    @property
    def meta(self):
        return self.api_schema.get_meta()

    @property
    def size(self):
        with open(
            os.path.join(self._model_path(self.model_id), MODEL_SIZE_FILE), "r"
        ) as f:
            return json.load(f)

    def example(self, n_samples, file_name=None, simple=True):
        eg = ExampleGenerator(model_id=self.model_id, config_json=self.config_json)
        return eg.example(n_samples=n_samples, file_name=file_name, simple=simple)
