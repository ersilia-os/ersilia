import os
import csv
import json
import time
import types
import tempfile
import importlib
import collections
import __main__ as main

from .. import logger
from ..serve.api import Api
from .session import Session
from datetime import datetime
from .base import ErsiliaBase
from ..lake.base import LakeBase
from ..utils import tmp_pid_file
from .modelbase import ModelBase
from ..serve.schema import ApiSchema
from ..utils.hdf5 import Hdf5DataLoader
from ..utils.csvfile import CsvDataLoader
from ..utils.terminal import yes_no_input
from ..serve.autoservice import AutoService
from ..io.output import TabularOutputStacker
from ..serve.standard_api import StandardCSVRunApi
from ..io.input import ExampleGenerator, BaseIOGetter
from .tracking import RunTracker, create_persistent_file
from ..io.readers.file import FileTyper, TabularFileReader
from ..utils.exceptions_utils.api_exceptions import ApiSpecifiedOutputError
from ..default import FETCHED_MODELS_FILENAME, MODEL_SIZE_FILE, CARD_FILE, EOS
from ..default import DEFAULT_BATCH_SIZE, APIS_LIST_FILE, INFORMATION_FILE
from ..utils.logging import make_temp_dir

try:
    import pandas as pd
except ModuleNotFoundError:
    pd = None


class ErsiliaModel(ErsiliaBase):
    def __init__(
        self,
        model,
        save_to_lake=True,
        service_class=None,
        config_json=None,
        credentials_json=None,
        verbose=None,
        fetch_if_not_available=True,
        preferred_port=None,
        track_runs=False,
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
        assert service_class in [
            None,
            "system",
            "venv",
            "conda",
            "docker",
            "pulled_docker",
            "hosted",
        ], "Wrong service class"
        self.service_class = service_class
        self.track_runs = track_runs
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
            try:
                do_fetch = yes_no_input(
                    "Requested model {0} is not available locally. Do you want to fetch it? [Y/n]".format(
                        self.model_id
                    ),
                    default_answer="Y",
                )
            except:
                self.logger.debug("Unable to capture user input. Fetching anyway.")
                do_fetch = True
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
        self.preferred_port = preferred_port
        self.autoservice = AutoService(
            model_id=self.model_id,
            service_class=self.service_class,
            config_json=self.config_json,
            preferred_port=preferred_port,
        )
        self._set_apis()
        self.session = Session(config_json=self.config_json)

        if track_runs:
            self._run_tracker = RunTracker(
                model_id=self.model_id, config_json=self.config_json
            )
        else:
            self._run_tracker = None

        self.logger.info("Done with initialization!")

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def is_valid(self):
        return self._is_valid

    def _set_api(self, api_name):
        # Don't want to override apis we explicitly write
        if hasattr(self, api_name):
            return

        def _method(input=None, output=None, batch_size=DEFAULT_BATCH_SIZE):
            return self.api(api_name, input, output, batch_size)

        setattr(self, api_name, _method)

    def _set_apis(self):
        apis_list = os.path.join(
            self._get_bundle_location(self.model_id), APIS_LIST_FILE
        )
        api_names = []
        if os.path.exists(apis_list):
            with open(apis_list, "r") as f:
                for l in f:
                    api_name = l.rstrip()
                    api_names += [api_name]
        if len(api_names) == 0:
            self.logger.debug("No apis found. Writing...")
            with open(apis_list, "w") as f:
                for (
                    api_name
                ) in self.autoservice.service._get_apis_from_where_available():
                    api_names += [api_name]
                    f.write(api_name + os.linesep)
        for api_name in api_names:
            self._set_api(api_name)
        self.apis_list = apis_list

    def _get_url(self):
        model_id = self.model_id
        tmp_file = tmp_pid_file(model_id)
        assert os.path.exists(
            tmp_file
        ), "Process ID file does not exist. Please serve the model first!"
        with open(tmp_file, "r") as f:
            for l in f:
                url = l.rstrip().split()[1]
        return url

    def _get_api_instance(self, api_name):
        url = self._get_url()
        if api_name is None:
            api_names = self.autoservice.get_apis()
            assert (
                len(api_names) == 1
            ), "More than one API found, please specificy api_name"
            api_name = api_names[0]
        api = Api(
            model_id=self.model_id,
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
            tmp_folder = make_temp_dir(prefix="ersilia-")
            is_h5_serializable = self.api_schema.is_h5_serializable(
                api_name=api.api_name
            )
            _temporary_prefix = "temporary"
            if is_h5_serializable:
                self.logger.debug("Output is HDF5 serializable")
                tmp_output = os.path.join(
                    tmp_folder, "{0}.h5".format(_temporary_prefix)
                )
            else:
                self.logger.debug("Output is not HDF5 serializable")
                tmp_output = os.path.join(
                    tmp_folder, "{0}.csv".format(_temporary_prefix)
                )
            for r in self._api_runner_iter(
                api=api, input=input, output=tmp_output, batch_size=batch_size
            ):
                continue
            if is_h5_serializable:
                data = Hdf5DataLoader()
            else:
                data = CsvDataLoader()
            data.load(tmp_output)
            if output == "numpy":
                return data.values[:]
            if output == "pandas":
                d = collections.OrderedDict()
                d["key"] = data.keys
                d["input"] = data.inputs
                df = pd.DataFrame(d)
                df[data.features] = data.values
                return df
            if output == "dict":
                d = collections.OrderedDict()
                d["keys"] = data.keys
                d["inputs"] = data.inputs
                d["values"] = data.values
                return d

    def _standard_api_runner(self, input, output):
        scra = StandardCSVRunApi(model_id=self.model_id, url=self._get_url())
        if not scra.is_ready():
            self.logger.debug(
                "Standard CSV Api runner is not ready for this particular model"
            )
            return None
        if not scra.is_amenable(input, output):
            self.logger.debug(
                "Standard CSV Api runner is not amenable for this model, input and output"
            )
            return None
        self.logger.debug("Starting standard runner")
        result = scra.post(input=input, output=output)
        return result

    @staticmethod
    def __output_is_file(output):
        if output is None:
            return False
        if type(output) != str:
            return False
        if output.endswith(".json"):
            return True
        if output.endswith(".csv"):
            return True
        if output.endswith(".tsv"):
            return True
        if output.endswith(".h5"):
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
        else:
            raise ApiSpecifiedOutputError
        if use_iter:
            return self._api_runner_iter
        else:
            return self._api_runner_return

    def _evaluate_do_cache_splits(self, input, output):
        if input is None:
            return False
        if output is None:
            return False
        if type(input) != str:
            return False
        if type(output) != str:
            return False
        if not os.path.exists(input):
            return False
        fti = FileTyper(input)
        if not fti.is_tabular():
            return False
        fto = FileTyper(output)
        if not fto.is_valid_output_file():
            return False
        return True

    def _do_cache_splits(self, input, output):
        self.tfr = None
        if self._evaluate_do_cache_splits(input, output):
            self.tfr = TabularFileReader(
                path=input,
                IO=BaseIOGetter(config_json=self.config_json).get(self.model_id),
            )
            if self.tfr.is_worth_splitting():
                return True
            else:
                self.tfr = None
                return False
        else:
            return False

    def api(
        self, api_name=None, input=None, output=None, batch_size=DEFAULT_BATCH_SIZE
    ):
        if self._do_cache_splits(input=input, output=output):
            splitted_inputs = self.tfr.split_in_cache()
            self.logger.debug("Split inputs:")
            self.logger.debug(" ".join(splitted_inputs))
            splitted_outputs = self.tfr.name_cached_output_files(
                splitted_inputs, output
            )
            for input_, output_ in zip(splitted_inputs, splitted_outputs):
                self.api_task(
                    api_name=api_name,
                    input=input_,
                    output=output_,
                    batch_size=batch_size,
                )
            TabularOutputStacker(splitted_outputs).stack(output)
            return output
        else:
            self.logger.debug("No file splitting necessary!")
            return self.api_task(
                api_name=api_name, input=input, output=output, batch_size=batch_size
            )

    def api_task(self, api_name, input, output, batch_size):
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

    def update_model_usage_time(self, model_id):
        file_name = os.path.join(EOS, FETCHED_MODELS_FILENAME)
        ts_str = str(time.time())
        with open(file_name, "r") as infile:
            models = dict(csv.reader(infile))
        infile.close()
        if model_id in models.keys():
            models[model_id] = ts_str

        with open(file_name, "w") as f:
            for key, values in models.items():
                f.write(f"{key},{values}\n")

    def serve(self):
        self.close()
        self.session.open(model_id=self.model_id, track_runs=self.track_runs)
        self.autoservice.serve()
        self.session.register_service_class(self.autoservice._service_class)
        self.url = self.autoservice.service.url
        self.pid = self.autoservice.service.pid
        self.scl = self.autoservice._service_class
        # self.update_model_usage_time(self.model_id) TODO: Check and reactivate

        # Start tracking to get the peak memory, memory usage and cpu time of the Model server(autoservice)
        if self._run_tracker is not None:
            create_persistent_file(self.model_id)
            memory_usage_serve, cpu_time_serve = self._run_tracker.get_memory_info()
            peak_memory_serve = self._run_tracker.get_peak_memory()

            session = Session(config_json=None)
            session.update_peak_memory(peak_memory_serve)
            session.update_total_memory(memory_usage_serve)
            session.update_cpu_time(cpu_time_serve)

    def close(self):
        self.autoservice.close()
        self.session.close()

    def get_apis(self):
        return self.autoservice.get_apis()

    def _run(
        self, input=None, output=None, batch_size=DEFAULT_BATCH_SIZE, track_run=False
    ):
        api_name = self.get_apis()[0]
        result = self.api(
            api_name=api_name, input=input, output=output, batch_size=batch_size
        )

        return result

    def _standard_run(self, input=None, output=None):
        t0 = time.time()
        t1 = None
        status_ok = False
        result = self._standard_api_runner(input=input, output=output)
        if type(output) is str:
            if os.path.exists(output):
                t1 = os.path.getctime(output)
        if t1 is not None:
            if t1 > t0:
                status_ok = True
        return result, status_ok

    def run(
        self,
        input=None,
        output=None,
        batch_size=DEFAULT_BATCH_SIZE,
        track_run=False,
        try_standard=True,
    ):
        self.logger.info("Starting runner")
        standard_status_ok = False
        if try_standard:
            self.logger.debug("Trying standard API")
            try:
                result, standard_status_ok = self._standard_run(
                    input=input, output=output
                )
            except Exception as e:
                self.logger.warning(
                    "Standard run did not work with exception {0}".format(e)
                )
                result = None
                standard_status_ok = False
                self.logger.debug("We will try conventional run.")
        if standard_status_ok:
            return result
        else:
            self.logger.debug("Trying conventional run")
            result = self._run(
                input=input, output=output, batch_size=batch_size, track_run=track_run
            )
        # Start tracking model run if track flag is used in serve
        if self._run_tracker is not None and track_run:
            self._run_tracker.track(input=input, result=result, meta=self._model_info)
            self._run_tracker.log(result=result, meta=self._model_info)
        return result

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

    def info(self):
        information_file = os.path.join(
            self._model_path(self.model_id), INFORMATION_FILE
        )
        with open(information_file, "r") as f:
            return json.load(f)

    @property
    def _model_info(self):
        return self.info()
