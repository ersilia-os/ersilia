import asyncio
import collections
import csv
import json
import os
import sys
import time
import types

from .. import logger
from ..default import (
    APIS_LIST_FILE,
    CARD_FILE,
    DEFAULT_BATCH_SIZE,
    EOS,
    FETCHED_MODELS_FILENAME,
    INFORMATION_FILE,
    MODEL_SIZE_FILE,
)
from ..hub.fetch.fetch import ModelFetcher
from ..io.input import BaseIOGetter, ExampleGenerator
from ..io.output import TabularOutputStacker
from ..io.readers.file import FileTyper, TabularFileReader
from ..serve.api import Api
from ..serve.autoservice import AutoService, PulledDockerImageService
from ..serve.schema import ApiSchema
from ..serve.standard_api import StandardCSVRunApi
from ..store.api import InferenceStoreApi
from ..store.utils import OutputSource
from ..utils import tmp_pid_file
from ..utils.csvfile import CsvDataLoader
from ..utils.exceptions_utils.api_exceptions import ApiSpecifiedOutputError
from ..utils.exceptions_utils.throw_ersilia_exception import throw_ersilia_exception
from ..utils.exceptions_utils.tracking_exceptions import TrackingNotSupportedError
from ..utils.hdf5 import Hdf5DataLoader
from ..utils.logging import make_temp_dir
from ..utils.terminal import yes_no_input
from .base import ErsiliaBase
from .modelbase import ModelBase
from .session import Session
from .tracking import RunTracker

try:
    import pandas as pd
except ModuleNotFoundError:
    pd = None


class ErsiliaModel(ErsiliaBase):
    """
    ErsiliaModel class for managing and interacting with different models.

    This class provides methods to fetch, serve, run, and close models form a model hub.
    It also supports tracking runs and handling various input and output formats.

    Parameters
    ----------
    model : str
        The identifier of the model.
    output_source : OutputSource, optional
        The source of the output, by default OutputSource.LOCAL_ONLY.
    service_class : str, optional
        The service class, by default None.
    config_json : dict, optional
        Configuration in JSON format, by default None.
    credentials_json : dict, optional
        Credentials in JSON format, by default None.
    verbose : bool, optional
        Verbosity flag, by default None.
    fetch_if_not_available : bool, optional
        Whether to fetch the model if not available locally, by default True.
    preferred_port : int, optional
        Preferred port for serving the model, by default None.
    track_runs : bool, optional
        Whether to track runs, by default False.
    cache: bool
        Whether to use redis cache or not
    maxmemory: float
        Fraction of memory used by redis

    Examples
    --------
    Fetching a model this requires to use asyncio since `fetch` is a coroutine.:

    .. code-block:: python

        model = ErsiliaModel(model="model_id")
        model.fetch()

    Serving a model:

    .. code-block:: python

        model = ErsiliaModel(model="model_id")
        model.serve()

    Running a model:

    .. code-block:: python

        model = ErsiliaModel(model="model_id")
        result = model.run(
            input="input_data.csv",
            output="output_data.csv",
        )

    Closing a model:

    .. code-block:: python

        model = ErsiliaModel(model="model_id")
        model.close()
    """

    def __init__(
        self,
        model: str,
        output_source: OutputSource = None,
        service_class: str = None,
        config_json: dict = None,
        credentials_json: dict = None,
        verbose: bool = None,
        fetch_if_not_available: bool = True,
        preferred_port: int = None,
        cache: bool = True,
        maxmemory: float = None,
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
            if hasattr(sys, "ps1"):
                self.logger.set_verbosity(0)

        assert service_class in [
            None,
            "system",
            "venv",
            "conda",
            "docker",
            "pulled_docker",
            "hosted",
        ], "Wrong service class"
        self.url = None
        self.pid = None
        self.service_class = service_class
        mdl = ModelBase(model)
        self._is_valid = mdl.is_valid()

        assert self._is_valid, (
            "The identifier {0} is not valid. Please visit the Ersilia Model Hub for valid identifiers".format(
                model
            )
        )
        self.config_json = config_json
        self.model_id = mdl.model_id
        self.slug = mdl.slug
        self.text = mdl.text
        self.output_source = output_source
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
                mf = ModelFetcher(
                    config_json=self.config_json, credentials_json=self.credentials_json
                )
                asyncio.run(mf.fetch(self.model_id))

        self.api_schema = ApiSchema(
            model_id=self.model_id, config_json=self.config_json
        )
        self.preferred_port = preferred_port
        self.autoservice = AutoService(
            model_id=self.model_id,
            service_class=self.service_class,
            config_json=self.config_json,
            preferred_port=preferred_port,
            cache=cache,
            maxmemory=maxmemory,
        )
        self._set_apis()
        self.logger.info("Done with initialization!")

    def fetch(self):
        """
        This method fetches the model from the Ersilia Model Hub.
        """
        mf = ModelFetcher(
            config_json=self.config_json, credentials_json=self.credentials_json
        )
        asyncio.run(mf.fetch(self.model_id))

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        This method is called when the runtime context is entered using the `with` statement.
        It starts serving the model.

        Returns
        -------
        ErsiliaModel
            The instance of the ErsiliaModel.
        """
        self.serve()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exit the runtime context related to this object.

        This method is called when the runtime context is exited using the `with` statement.
        It closes the model services and session.

        Parameters
        ----------
        exc_type : type
            The exception type.
        exc_val : Exception
            The exception instance.
        exc_tb : traceback
            The traceback object.
        """
        self.close()

    def is_valid(self):
        """
        Check if the model identifier is valid.

        This method verifies if the provided model identifier is valid by checking its existence
        and validity in the model hub.

        Returns
        -------
        bool
            True if the model identifier is valid, False otherwise.
        """
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
        assert os.path.exists(tmp_file), (
            "Process ID file does not exist. Please serve the model first!"
        )
        with open(tmp_file, "r") as f:
            for l in f:
                url = l.rstrip().split()[1]
        return url

    def _get_api_instance(self, api_name):
        url = self._get_url()
        if api_name is None:
            api_names = self.autoservice.get_apis()
            assert len(api_names) == 1, (
                "More than one API found, please specificy api_name"
            )
            api_name = api_names[0]
        api = Api(
            model_id=self.model_id,
            url=url,
            api_name=api_name,
            config_json=self.config_json,
        )
        return api

    def _api_runner_iter(self, api: Api, input: str, output: str, batch_size: int):
        """
        Run the API in an iterative manner.

        This method executes the API in an iterative manner, yielding results one by one.
        It is useful for handling large datasets that need to be processed in batches.

        Parameters
        ----------
        api : Api
            The API instance to run.
        input : str
            The input data.
        output : str
            The output data.
        batch_size : int
            The batch size.

        Yields
        ------
        Any
            The result of each API call.
        """
        for result in api.post(input=input, output=output, batch_size=batch_size):
            assert result is not None, (
                "Something went wrong. Please contact us at hello@ersila.io"
            )
            yield result

    def _api_runner_return(self, api: Api, input: str, output: str, batch_size: int):
        """
        Run the API and return the results.

        This method executes the API and returns the results in the specified output format.
        It handles different output formats such as JSON, HDF5, CSV, NumPy arrays, Pandas DataFrames, and dictionaries.

        Parameters
        ----------
        api : Api
            The API instance to run.
        input : str
            The input data.
        output : str
            The output data.
        batch_size : int
            The batch size.

        Returns
        -------
        Any
            The result of the API run in the specified output format.
        """
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

    def _standard_api_runner(self, input, output, batch_size):
        scra = StandardCSVRunApi(model_id=self.model_id, url=self._get_url())
        if not scra.is_amenable(output):
            self.logger.debug(
                "Standard CSV Api runner is not amenable for this model, input and output"
            )
            return None
        self.logger.debug("Starting standard runner")
        result = scra.post(
            input=input,
            output=output,
            batch_size=batch_size,
            output_source=self.output_source,
        )
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
        """
        Run the specified API with the given input and output.

        This method executes the specified API(usually with the end point run) using the provided input and output parameters.
        It handles file splitting and caching if necessary.

        Parameters
        ----------
        api_name : str, optional
            The name of the API to run, by default None.
        input : str, optional
            The input data, by default None.
        output : str, optional
            The output data, by default None.
        batch_size : int, optional
            The batch size, by default DEFAULT_BATCH_SIZE.

        Returns
        -------
        Any
            The result of the API run.
        """
        if OutputSource.is_precalculation_enabled(self.output_source):
            store = InferenceStoreApi(model_id=self.model_id, output=output)
            return store.get_precalculations(input)
        elif self._do_cache_splits(input=input, output=output):
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
        """
        Run the specified API task with the given input and output.

        This method executes the specified API task using the provided input and output parameters.
        It returns the result of the API task, which can be a generator, file, or other data types.

        Parameters
        ----------
        api_name : str
            The name of the API to run.
        input : str
            The input data.
        output : str
            The output data.
        batch_size : int
            The batch size.

        Returns
        -------
        Any
            The result of the API task.
        """
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
        """
        Update the model usage time.

        This method updates the usage time of the specified model by recording the current timestamp
        in the fetched models file.

        Parameters
        ----------
        model_id : str
            The identifier of the model.
        """
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

    def setup(self):
        """
        Setup the necessary requirements for the model.

        This method ensures that the required dependencies and resources for the model are available.
        """
        pass  # TODO Implement whenever this is necessary

    @throw_ersilia_exception()
    def serve(self, track_runs=None):
        """
        Serve the model by starting the necessary services.

        This method sets up the required dependencies, opens a session, and starts the model service.
        It registers the service class and output source, updates the model's URL and process ID (PID),
        and tracks resource usage if tracking is enabled.

        Parameters
        ----------
        track_runs : str, optional
            Whether to track runs, by default is None (i.e. do not track).
        """
        self.logger.debug("Starting serve")
        self.session = Session(config_json=self.config_json)
        self.run_tracker = None
        self.tracking_use_case = track_runs
        self.track = False
        if track_runs is not None:
            if not isinstance(self.autoservice.service, PulledDockerImageService):
                self.logger.warning(
                    "Tracking runs is currently only supported for Dockerized models"
                )
                raise TrackingNotSupportedError()
            else:
                self.track = True
                self.run_tracker = RunTracker(
                    model_id=self.model_id,
                    config_json=self.config_json,
                    use_case=track_runs,
                )
        self.setup()
        self.close()
        self.session.open(model_id=self.model_id, track_runs=self.track)
        self.autoservice.serve()
        self.session.register_service_class(self.autoservice._service_class)
        self.session.register_output_source(self.output_source)
        if self.track:
            self.session.register_tracking_use_case(self.tracking_use_case)
        self.url = self.autoservice.service.url
        self.pid = self.autoservice.service.pid
        self.scl = self.autoservice._service_class
        self.logger.debug("Done with basic session registration")

    def close(self):
        """
        Close the model services and session.

        This method stops the model service and closes the session.
        """
        self.session = Session(config_json=self.config_json)
        self.logger.debug("Closing session {0}".format(self.session._session_dir))
        self.logger.debug("Stopping service")
        self.autoservice.close()
        self.session.close()

    def get_apis(self):
        """
        Get the list of available APIs for the model.

        This method retrieves the list of APIs that are available for the model.

        Returns
        -------
        list
            The list of available APIs.
        """
        return self.autoservice.get_apis()

    def _run(self, input=None, output=None, batch_size=DEFAULT_BATCH_SIZE):
        api_name = self.get_apis()[0]
        result = self.api(
            api_name=api_name, input=input, output=output, batch_size=batch_size
        )

        return result

    def _standard_run(self, input=None, output=None, batch_size=DEFAULT_BATCH_SIZE):
        t0 = time.time()
        t1 = None
        status_ok = False
        result = self._standard_api_runner(
            input=input, output=output, batch_size=batch_size
        )
        if type(output) is str:  # TODO Redundant, should be removed
            if os.path.exists(output):
                t1 = os.path.getctime(output)
        if t1 is not None:
            if t1 > t0:  # TODO Why would this be less?
                status_ok = True
        return (
            result,
            status_ok,
        )  # TODO This doesn't actually check whether the result exists and isn't null

    @throw_ersilia_exception()
    def run(
        self,
        input: str = None,
        output: str = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        """
        Run the model with the given input and output.

        This method executes the model using the provided input and output parameters.
        It first tries to run the model using the standard API, and if that fails, it falls back to
        the conventional run method. It also tracks the run if the track_run flag is set.

        Parameters
        ----------
        input : str, optional
            The input data, by default None.
        output : str, optional
            The output data, by default None.
        batch_size : int, optional
            The batch size, by default DEFAULT_BATCH_SIZE.

        Returns
        -------
        Any
            The result of the model run(such as output csv file name, json).
        """
        t0 = time.time()

        session = Session(config_json=self.config_json)
        track_run = session.tracking_status()

        # Currently, tracking is only accepted for CSV files
        if output is None:
            if track_run:
                raise TrackingNotSupportedError()
        else:
            if not output.endswith(".csv") and track_run:
                raise TrackingNotSupportedError()

        # Init the run tracker
        if track_run:
            use_case = session.get_tracking_use_case()
            self.logger.debug(
                "Initializing the run trackers (use case {0})".format(use_case)
            )
            self.run_tracker = RunTracker(
                model_id=self.model_id, config_json=self.config_json, use_case=use_case
            )
        self.logger.info("Starting runner")

        # TODO The logic should be in a try except else finally block
        standard_status_ok = False
        self.logger.debug("Trying standard API")
        try:
            result, standard_status_ok = self._standard_run(
                input=input, output=output, batch_size=batch_size
            )
        except Exception as e:
            self.logger.warning(
                "Standard run did not work with exception {0}".format(e)
            )
            result = None
            standard_status_ok = False
            self.logger.debug("We will try conventional run.")

        if not standard_status_ok:
            self.logger.debug("Trying conventional run")
            if track_run:
                self.logger.error(
                    "With conventional runner tracker will not be enabled. Disable tracking at serve time if you want to proceed."
                )
                raise TrackingNotSupportedError()
            result = self._run(input=input, output=output, batch_size=batch_size)

        # Collect metrics sampled during run if tracking is enabled
        if track_run:
            self.logger.debug("Collecting metrics")
            model_info = self.info()
            if "metadata" in model_info:
                metadata = model_info["metadata"]
            elif "card" in model_info:
                metadata = model_info["card"]
            else:
                metadata = {}

            t1 = time.time()
            time_taken = t1 - t0
            self.run_tracker.track(input, output, metadata, time_taken)

        return result

    @property
    def paths(self):
        """
        Get the paths related to the model.

        This property returns a dictionary containing various paths related to the model,
        such as the destination path, repository path, and BentoML path.

        Returns
        -------
        dict
            The dictionary containing paths.
        """
        p = {
            "dest": self._model_path(self.model_id),
            "repository": self._get_bundle_location(self.model_id),
            "bentoml": self._get_bentoml_location(self.model_id),
        }
        return p

    @property
    def input_type(self):
        """
        Get the input type of the model.

        This property reads the input type information from the model's card file and returns it
        as a list of input types.

        Returns
        -------
        list
            The list of input types(such as compounds).
        """
        with open(os.path.join(self._model_path(self.model_id), CARD_FILE), "r") as f:
            return [x.lower() for x in json.load(f)["Input"]]

    @property
    def output_type(self):
        """
        Get the output type of the model.

        This property reads the output type information from the model's card file and returns it
        as a list of output types.

        Returns
        -------
        list
            The list of output types(such as Descriptor, score, probability etc...).
        """
        with open(os.path.join(self._model_path(self.model_id), CARD_FILE), "r") as f:
            return [x.lower() for x in json.load(f)["Output"]]

    @property
    def schema(self):
        """
        Get the schema of the model.

        This property returns the schema of the model, which defines the structure and format
        of the model's input and output data.

        Returns
        -------
        dict
            The schema of the model.
        """
        return self.api_schema.schema

    @property
    def meta(self):
        """
        Get the metadata of the model.

        This property returns the metadata of the model, which provides additional information
        about the model, such as its description, version, and author.

        Returns
        -------
        dict
            The metadata of the model.
        """
        return self.api_schema.get_meta()

    @property
    def size(self):
        """
        Get the size of the model.

        This property reads the size information from the model's size file and returns it
        as a dictionary.

        Returns
        -------
        dict
            The size of the model.
        """
        with open(
            os.path.join(self._model_path(self.model_id), MODEL_SIZE_FILE), "r"
        ) as f:
            return json.load(f)

    def example(self, n_samples, file_name=None, simple=True):
        """
        Generate example data for the model.

        This method generates example data for the model using the specified number of samples.
        The generated data can be saved to a file if a file name is provided.

        Parameters
        ----------
        n_samples : int
            The number of samples to generate.
        file_name : str, optional
            The file name to save the examples, by default None.
        simple : bool, optional
            Whether to generate simple examples, by default True.

        Returns
        -------
        Any
            The generated example data(path, list of smiles etc...).
        """
        eg = ExampleGenerator(model_id=self.model_id, config_json=self.config_json)
        return eg.example(n_samples=n_samples, file_name=file_name, simple=simple)

    def info(self):
        """
        Get the information of the model.

        This method reads the information file of the model and returns its content as a dictionary.

        Returns
        -------
        dict
            The information of the model.
        """
        information_file = os.path.join(
            self._model_path(self.model_id), INFORMATION_FILE
        )
        with open(information_file, "r") as f:
            return json.load(f)
