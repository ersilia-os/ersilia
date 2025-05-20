import importlib
import json
import os
import time
import uuid

import docker
import requests

from .. import ErsiliaBase, throw_ersilia_exception
from ..db.environments.localdb import EnvironmentDb
from ..db.environments.managers import DockerManager
from ..default import (
    ALLOWED_API_NAMES,
    APIS_LIST_FILE,
    DEFAULT_DOCKER_NETWORK_BRIDGE,
    DEFAULT_DOCKER_NETWORK_NAME,
    DEFAULT_VENV,
    DOCKERHUB_ORG,
    INFORMATION_FILE,
    IS_FETCHED_FROM_HOSTED_FILE,
    PACK_METHOD_BENTOML,
    PACK_METHOD_FASTAPI,
    PACKMODE_FILE,
)
from ..setup.requirements.bentoml_requirement import BentoMLRequirement

# from ..setup.requirements.bentoml import BentoMLRequirement
from ..setup.requirements.conda import CondaRequirement
from ..setup.requirements.docker import DockerRequirement
from ..tools.bentoml.exceptions import BentoMLException
from ..utils.conda import SimpleConda, StandaloneConda
from ..utils.docker import SimpleDocker, model_image_version_reader, set_docker_host
from ..utils.exceptions_utils.serve_exceptions import (
    BadGatewayError,
    DockerNotActiveError,
)
from ..utils.logging import make_temp_dir
from ..utils.ports import find_free_port
from ..utils.terminal import run_command
from ..utils.venv import SimpleVenv

SLEEP_SECONDS = 1
TIMEOUT_SECONDS = 1000


class BaseServing(ErsiliaBase):
    """
    Base class for serving models.

    Parameters
    ----------
    model_id : str
        The ID of the model to be served.
    config_json : dict, optional
        Configuration settings in JSON format.
    preferred_port : int, optional
        Preferred port for serving the model.
    url : str, optional
        URL for the served model.
    """

    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.bundle_tag = self._get_latest_bundle_tag(model_id=self.model_id)
        self.port = preferred_port

    def _get_apis_from_apis_list(self):
        self.logger.debug("Getting APIs from list file")
        file_name = os.path.join(
            self._get_bundle_location(self.model_id), APIS_LIST_FILE
        )
        if not os.path.exists(file_name):
            return None
        with open(file_name, "r") as f:
            apis_list = []
            for l in f:
                apis_list += [l.rstrip()]
        if len(apis_list) > 0:
            return apis_list
        else:
            return None

    def _get_info_from_bento(self):
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "information.json")
        cmd = "bentoml info --quiet {0}:{1} > {2}".format(
            self.model_id, self.bundle_tag, tmp_file
        )
        self.logger.debug(
            "Getting info from BentoML and storing in {0}".format(tmp_file)
        )
        # Check command success first
        result = run_command(cmd)
        if result.returncode != 0:
            raise BentoMLException(f"BentoML info failed: {result.stderr}")

        # Handle JSON parsing errors here
        try:
            with open(tmp_file, "r") as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            self.logger.error(f"Invalid BentoML output: {e}")
            raise BentoMLException("Corrupted BentoML installation detected") from e

    def _get_apis_from_bento(self):
        self.logger.debug("Getting APIs from Bento")
        bento_requirement = BentoMLRequirement()

        try:
            info = self._get_info_from_bento()
        except BentoMLException as e:
            # Handle both command failures and JSON errors here
            if "Corrupted BentoML installation" in str(e):
                self.logger.warning("Attempting BentoML cleanup...")
                try:
                    bento_requirement._cleanup_corrupted_bentoml()
                    self.logger.info("Retrying API fetch after cleanup")
                    info = self._get_info_from_bento()  # Retry
                except Exception as cleanup_error:
                    raise BentoMLException(
                        f"Cleanup failed: {cleanup_error}"
                    ) from cleanup_error
            else:
                raise  # Re-raise unrelated errors

        try:
            return [item["name"] for item in info["apis"]]
        except KeyError as e:
            raise BentoMLException(f"Invalid API format: {e}") from e

    def _get_apis_from_fastapi(self):
        bundle_path = self._model_path(self.model_id)
        apis_list = []
        framework_path = os.path.join(bundle_path, "model", "framework")
        if not os.path.exists(framework_path):
            return None
        for fn in os.listdir(framework_path):
            if fn.endswith(".sh"):
                api_name = fn.split(".")[0]
                apis_list += [api_name]
        return apis_list

    def _get_apis_from_where_available(self):
        apis_list = self._get_apis_from_apis_list()
        if apis_list is None:
            pack_method = self._resolve_pack_method_source(self.model_id)
            if pack_method == PACK_METHOD_FASTAPI:
                self.logger.debug("Getting APIs from FastAPI")
                apis_list = self._get_apis_from_fastapi()
            elif pack_method == PACK_METHOD_BENTOML:
                self.logger.debug("Getting APIs from BentoML")
                apis_list = self._get_apis_from_bento()
            else:
                raise
        if apis_list is None:
            apis_list = []
        for api in apis_list:
            yield api

    def _api_with_url(self, api_name, input):
        """
        Call an API with the given URL and input.

        Parameters
        ----------
        api_name : str
            Name of the API to call.
        input : dict
            Input data for the API.

        Returns
        -------
        dict
            Response from the API.
        """
        if self.url is None:
            return
        self.logger.debug("Using URL: {0}".format(self.url))
        response = requests.post("{0}/{1}".format(self.url, api_name), json=input)
        return response.json()


class _BentoMLService(BaseServing):
    def __init__(self, model_id, config_json=None, preferred_port=None):
        BaseServing.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )
        self.SEARCH_PRE_STRING = "* Running on "
        self.SEARCH_SUF_STRING = "Press CTRL+C to quit"
        self.ERROR_STRING = "error"

    def serve(self, runcommand_func=None):
        """
        Serve the model using BentoML.

        Parameters
        ----------
        runcommand_func : function, optional
            Function to run the command. If None, the command is run from the shell.
        """
        self.logger.debug("Trying to serve model with BentoML locally")
        preferred_port = self.port
        self.port = find_free_port(preferred_port=preferred_port)
        if self.port != preferred_port:
            self.logger.warning(
                "Port {0} was already in use. Using {1} instead".format(
                    preferred_port, self.port
                )
            )
        self.logger.debug("Free port: {0}".format(self.port))
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_script = os.path.join(tmp_folder, "serve.sh")
        tmp_file = os.path.join(tmp_folder, "serve.log")
        tmp_pid = os.path.join(tmp_folder, "serve.pid")
        sl = ["#!/bin/bash"]
        sl += [
            "bentoml serve {0}:{1} --port {2} &> {3} &".format(
                self.model_id, self.bundle_tag, self.port, tmp_file
            )
        ]
        sl += ["_pid=$!"]
        sl += ['echo "$_pid" > {0}'.format(tmp_pid)]
        self.logger.debug("Writing on {0}".format(tmp_script))
        with open(tmp_script, "w") as f:
            for l in sl:
                self.logger.debug(l)
                f.write(l + os.linesep)
        cmd = "bash {0}".format(tmp_script)
        if runcommand_func is None:
            self.logger.debug("Run command function not available. Running from shell")
            run_command(cmd)
        else:
            self.logger.debug("Run command function available")
            runcommand_func(cmd)
        with open(tmp_pid, "r") as f:
            self.pid = int(f.read().strip())
            self.logger.debug("Process id: {0}".format(self.pid))
        _logged_file_done = False
        _logged_server_done = False
        for it in range(int(TIMEOUT_SECONDS / SLEEP_SECONDS)):
            self.logger.debug("Trying to wake up. Iteration: {0}".format(it))
            self.logger.debug(
                "Timeout: {0} Sleep time: {1}".format(TIMEOUT_SECONDS, SLEEP_SECONDS)
            )
            if not os.path.exists(tmp_file):
                if not _logged_file_done:
                    self.logger.debug("Waiting for file {0}".format(tmp_file))
                _logged_file_done = True
                time.sleep(SLEEP_SECONDS)
                continue
            self.logger.debug("Temporary file available: {0}".format(tmp_file))
            # If error string is identified, finish
            with open(tmp_file, "r") as f:
                r = f.read()
                if self.ERROR_STRING in r.lower():
                    self.logger.warning("Error string found in: {0}".format(r))
                    # TODO perhaps find a better error string.
                    # self.url = None
                    # return
            self.logger.debug("No error strings found in temporary file")
            # If everything looks good, wait until server is ready
            with open(tmp_file, "r") as f:
                r = f.read()
                if self.SEARCH_PRE_STRING not in r or self.SEARCH_SUF_STRING not in r:
                    if not _logged_server_done:
                        self.logger.debug("Waiting for server")
                    else:
                        self.logger.debug("Server logging done")
                    time.sleep(SLEEP_SECONDS)
                    _logged_server_done = True
                    continue
            self.logger.debug("Server is ready. Trying to get URL")
            # When the search strings are found get url
            with open(tmp_file, "r") as f:
                for l in f:
                    if self.SEARCH_PRE_STRING in l:
                        self.url = (
                            l.split(self.SEARCH_PRE_STRING)[1].split(" ")[0].rstrip()
                        )
                        self.logger.debug("URL found: {0}".format(self.url))
                        return
                self.logger.debug("Search strings not found yet")
        self.logger.debug("No URL found")
        self.url = None

    def close(self):
        """
        Close the BentoML service by killing the process.
        """
        try:
            os.kill(self.pid, 9)
        except:
            self.logger.info("PID {0} is unassigned".format(self.pid))


class _FastApiService(BaseServing):
    def __init__(self, model_id, config_json=None, preferred_port=None):
        BaseServing.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )
        self.SEARCH_PRE_STRING = "Uvicorn running on "
        self.SEARCH_SUF_STRING = "(Press CTRL+C to quit)"
        self.ERROR_STRING = "error"
        self.conda = SimpleConda()

    def serve(self, runcommand_func=None):
        """
        Serve the model using FastAPI.

        Parameters
        ----------
        runcommand_func : function, optional
            Function to run the command. If None, the command is run from the shell.
        """
        bundle_path = self._get_bundle_location(self.model_id)
        self.logger.debug("Trying to serve model with FastAPI locally")
        preferred_port = self.port
        self.port = find_free_port(preferred_port=preferred_port)
        if self.port != preferred_port:
            self.logger.warning(
                "Port {0} was already in use. Using {1} instead".format(
                    preferred_port, self.port
                )
            )
        self.logger.debug("Free port: {0}".format(self.port))
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_script = os.path.join(tmp_folder, "serve.sh")
        tmp_file = os.path.join(tmp_folder, "serve.log")
        tmp_pid = os.path.join(tmp_folder, "serve.pid")
        sl = [
            "ersilia_model_serve --bundle_path {0} --port {1} &> {2} &".format(
                bundle_path, self.port, tmp_file
            )
        ]
        sl += ["_pid=$!"]
        sl += ['echo "$_pid" > {0}'.format(tmp_pid)]
        self.logger.debug("Writing on {0}".format(tmp_script))
        self.conda.create_executable_bash_script(
            environment=self.model_id, commandlines=sl, file_name=tmp_script
        )
        cmd = "bash {0}".format(tmp_script)
        if runcommand_func is None:
            self.logger.debug("Run command function not available. Running from shell")
            run_command(cmd)
        else:
            self.logger.debug("Run command function available")
            runcommand_func(cmd)
        with open(tmp_pid, "r") as f:
            self.pid = int(f.read().strip())
            self.logger.debug("Process id: {0}".format(self.pid))
        _logged_file_done = False
        _logged_server_done = False
        for it in range(int(TIMEOUT_SECONDS / SLEEP_SECONDS)):
            self.logger.debug("Trying to wake up. Iteration: {0}".format(it))
            self.logger.debug(
                "Timeout: {0} Sleep time: {1}".format(TIMEOUT_SECONDS, SLEEP_SECONDS)
            )
            if not os.path.exists(tmp_file):
                if not _logged_file_done:
                    self.logger.debug("Waiting for file {0}".format(tmp_file))
                _logged_file_done = True
                time.sleep(SLEEP_SECONDS)
                continue
            self.logger.debug("Temporary file available: {0}".format(tmp_file))
            # If error string is identified, finish
            with open(tmp_file, "r") as f:
                r = f.read()
                if self.ERROR_STRING in r.lower():
                    self.logger.warning("Error string found in: {0}".format(r))
                    # TODO perhaps find a better error string.
                    # self.url = None
                    # return
            self.logger.debug("No error strings found in temporary file")
            # If everything looks good, wait until server is ready
            with open(tmp_file, "r") as f:
                r = f.read()
                if self.SEARCH_PRE_STRING not in r or self.SEARCH_SUF_STRING not in r:
                    if not _logged_server_done:
                        self.logger.debug("Waiting for server")
                    else:
                        self.logger.debug("Server logging done")
                    time.sleep(SLEEP_SECONDS)
                    _logged_server_done = True
                    continue
            self.logger.debug("Server is ready. Trying to get URL")
            # When the search strings are found get url
            with open(tmp_file, "r") as f:
                for l in f:
                    if self.SEARCH_PRE_STRING in l:
                        self.url = (
                            l.split(self.SEARCH_PRE_STRING)[1].split(" ")[0].rstrip()
                        )
                        self.logger.debug("URL found: {0}".format(self.url))
                        return
                self.logger.debug("Search strings not found yet")
        self.logger.debug("No URL found")
        self.url = None

    def close(self):
        """
        Close the FastAPI service by killing the process.
        """
        try:
            os.kill(self.pid, 9)
        except:
            self.logger.info("PID {0} is unassigned".format(self.pid))


class _LocalService(ErsiliaBase):
    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        self.model_id = model_id
        ErsiliaBase.__init__(self, config_json=config_json)
        pack_method = self._resolve_pack_method_source(self.model_id)
        self.logger.debug("Pack method is: {0}".format(pack_method))
        if pack_method == PACK_METHOD_FASTAPI:
            self.server = _FastApiService(
                model_id,
                config_json=config_json,
                preferred_port=preferred_port,
            )
        elif pack_method == PACK_METHOD_BENTOML:
            self.server = _BentoMLService(
                model_id,
                config_json=config_json,
                preferred_port=preferred_port,
            )
        else:
            raise Exception("Model is not a valid BentoML or FastAPI model")

    def _get_apis_from_where_available(self):
        return self.server._get_apis_from_where_available()

    def local_serve(self, runcommand_func=None):
        """
        Serve the model locally.

        Parameters
        ----------
        runcommand_func : function, optional
            Function to run the command. If None, the command is run from the shell.
        """
        self.server.serve(runcommand_func=runcommand_func)
        self.url = self.server.url
        self.pid = self.server.pid
        self.port = self.server.port

    def local_close(self):
        """
        Close the local service.
        """
        self.server.close()


class SystemBundleService(_LocalService):
    """
    Service class for managing system bundle models.

    Parameters
    ----------
    model_id : str
        The ID of the model to be served.
    config_json : dict, optional
        Configuration settings in JSON format.
    preferred_port : int, optional
        Preferred port for serving the model.
    url : str, optional
        URL for the served model.
    """

    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        _LocalService.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        Returns
        -------
        SystemBundleService
            The instance of the service.
        """
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Exit the runtime context related to this object.

        Parameters
        ----------
        exception_type : type
            The exception type.
        exception_value : Exception
            The exception instance.
        traceback : traceback
        """
        self.close()

    def _run_command(self, cmd):
        return run_command(cmd)

    def is_available(self):
        """
        Check if the system bundle service is available.

        Returns
        -------
        bool
            True if the service is available, False otherwise.
        """
        fn = os.path.join(self._dest_dir, self.model_id, PACKMODE_FILE)
        avail = False
        if os.path.exists(fn):
            with open(fn, "r") as f:
                pack_mode = f.read()
                if pack_mode == "system":
                    avail = True
        self.url = None
        return avail

    def serve(self):
        """
        Serve the model using the system bundle service.
        """
        self.local_serve()

    def close(self):
        """
        Close the system bundle service.
        """
        self.local_close()

    def api(self, api_name, input):
        """
        Call an API with the given name and input.

        Parameters
        ----------
        api_name : str
            Name of the API to call.
        input : dict
            Input data for the API.

        Returns
        -------
        dict
            Response from the API.
        """
        return self._api_with_url(api_name, input)


class VenvEnvironmentService(_LocalService):
    """
    Service class for managing virtual environment models.

    Parameters
    ----------
    model_id : str
        The ID of the model to be served.
    config_json : dict, optional
        Configuration settings in JSON format.
    preferred_port : int, optional
        Preferred port for serving the model.
    url : str, optional
        URL for the served model.
    """

    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        _LocalService.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )
        self.venv = SimpleVenv(self._model_path(model_id))

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        Returns
        -------
        VenvEnvironmentService
            The instance of the service.
        """
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Exit the runtime context related to this object.

        Parameters
        ----------
        exception_type : type
            The exception type.
        exception_value : Exception
            The exception instance.
        traceback : traceback
            The traceback object.
        """
        self.close()

    def _model_path(self, model_id):
        return os.path.join(self._dest_dir, model_id)

    def _run_command(self, cmd):
        return self.venv.run_commandlines(DEFAULT_VENV, cmd)

    def is_available(self):
        """
        Check if the virtual environment service is available.

        Returns
        -------
        bool
            True if the service is available, False otherwise.
        """
        if not self.venv.exists(DEFAULT_VENV):
            return False

    def serve(self):
        """
        Serve the model using the virtual environment service.
        """
        self.local_serve(self._run_command)

    def close(self):
        """
        Close the virtual environment service.
        """
        self.local_close()

    def api(self, api_name, input):
        """
        Call an API with the given name and input.

        Parameters
        ----------
        api_name : str
            Name of the API to call.
        input : dict
            Input data for the API.

        Returns
        -------
        dict
            Response from the API.
        """
        return self._api_with_url(api_name, input)


class CondaEnvironmentService(_LocalService):
    """
    Service class for managing Conda environment models.

    Parameters
    ----------
    model_id : str
        The ID of the model to be served.
    config_json : dict, optional
        Configuration settings in JSON format.
    preferred_port : int, optional
        Preferred port for serving the model.
    url : str, optional
        URL for the served model.
    """

    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        _LocalService.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )
        if CondaEnvironmentService.is_single_model_without_conda():
            self._is_singe = True
            self.conda = StandaloneConda()
        else:
            self._is_singe = False
            self.db = EnvironmentDb()
            self.db.table = "conda"
            self.conda = SimpleConda()

    @staticmethod
    def is_single_model_without_conda():
        """
        Check if there is a single model without conda.

        Returns
        -------
        bool
            True if conda is not installed, False otherwise.
        """
        conda_checker = CondaRequirement()
        # Returns True if conda is not installed and False otherwise
        return not conda_checker.is_installed()

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        Returns
        -------
        CondaEnvironmentService
            The instance of the service.
        """
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Exit the runtime context related to this object.

        Parameters
        ----------
        exception_type : type
            The exception type.
        exception_value : Exception
            The exception instance.
        traceback : traceback
        """
        self.close()

    def _get_env_name(self):
        if self._is_singe:
            return self.model_id
        else:
            envs = list(self.db.envs_of_model(self.model_id))
            for env in envs:
                if self.conda.exists(env):
                    return env
            return None

    def _run_command(self, cmd):
        env = self._get_env_name()
        self.logger.debug("Running on conda environment {0}".format(env))
        return self.conda.run_commandlines(env, cmd)

    def is_available(self):
        """
        Check if the Conda environment service is available.

        Returns
        -------
        bool
            True if the service is available, False otherwise.
        """
        env = self._get_env_name()
        if env is not None:
            return True
        else:
            return False

    def serve(self):
        """
        Serve the model using the Conda environment service.
        """
        self.local_serve(self._run_command)

    def close(self):
        """
        Close the Conda environment service.
        """
        self.local_close()

    def api(self, api_name, input):
        """
        Call an API with the given name and input.

        Parameters
        ----------
        api_name : str
            Name of the API to call.
        input : dict
            Input data for the API.

        Returns
        -------
        dict
            Response from the API.
        """
        return self._api_with_url(api_name, input)


class DockerImageService(BaseServing):
    """
    Service class for managing Docker image models.

    Parameters
    ----------
    model_id : str
        The ID of the model to be served.
    config_json : dict, optional
        Configuration settings in JSON format.
    preferred_port : int, optional
        Preferred port for serving the model.
    url : str, optional
        URL for the served model.
    """

    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        self._is_docker_active()
        BaseServing.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )
        self.db = EnvironmentDb()
        self.db.table = "docker"
        self.docker = SimpleDocker()
        self.dm = DockerManager(config_json=config_json, preferred_port=preferred_port)
        self.pid = -1

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        Returns
        -------
        DockerImageService
            The instance of the service.
        """
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Exit the runtime context related to this object.

        Parameters
        ----------
        exception_type : type
            The exception type.
        exception_value : Exception
            The exception instance.
        traceback : traceback
        """
        self.close()

    def _get_env_name(self):
        envs = list(self.db.envs_of_model(self.model_id))
        for env in envs:
            org, img, tag = self.docker._splitter(env)
            if self.docker.exists(org=org, img=img, tag=tag):
                self.logger.debug("Docker image found {0}".format(env))
                return env
        return None

    @throw_ersilia_exception()
    def _is_docker_active(self):
        dr = DockerRequirement()
        if not dr.is_active():
            raise DockerNotActiveError()
        return True

    def is_available(self):
        """
        Check if the Docker image service is available.

        Returns
        -------
        bool
            True if the service is available, False otherwise.
        """
        env = self._get_env_name()
        if env is not None:
            self.logger.debug("Docker image service available")
            return True
        else:
            self.logger.debug("Docker image service not available")
            return False

    def serve(self):
        """
        Serve the model using the Docker image service.
        """
        self.logger.debug("Calling docker manager")
        res = self.dm.run(self.model_id)
        self.container_name = res["container_name"]
        self.port = res["port"]
        self.url = "http://0.0.0.0:{0}".format(self.port)

    def close(self):
        """
        Close the Docker image service.
        """
        self.df.stop_containers(self.model_id)

    def api(self, api_name, input):
        """
        Call an API with the given name and input.

        Parameters
        ----------
        api_name : str
            Name of the API to call.
        input : dict
            Input data for the API.

        Returns
        -------
        dict
            Response from the API.
        """
        return self._api_with_url(api_name, input)


# TODO: Include 'pip' within available service_class
class PipInstalledService(BaseServing):
    """
    Service class for managing pip-installed models.

    Parameters
    ----------
    model_id : str
        The ID of the model to be served.
    config_json : dict, optional
        Configuration settings in JSON format.
    preferred_port : int, optional
        Preferred port for serving the model.
    url : str, optional
        URL for the served model.
    """

    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        BaseServing.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )
        self.pid = -1

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        Returns
        -------
        PipInstalledService
            The instance of the service.
        """
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Exit the runtime context related to this object.

        Parameters
        ----------
        exception_type : type
            The exception type.
        exception_value : Exception
            The exception instance.
        traceback : traceback
        """
        self.close()

    def _import(self):
        try:
            model = importlib.import_module(self.model_id, package=None)
            return model
        except ModuleNotFoundError:
            return None

    def is_available(self):
        """
        Check if the pip-installed service is available.

        Returns
        -------
        bool
            True if the service is available, False otherwise.
        """
        model = self._import()
        if model is not None:
            return True
        else:
            return False

    def serve(self):
        """
        Serve the model using the pip-installed service.
        """
        model = self._import()
        self.mdl = model.load()

    def close(self):
        """
        Close the pip-installed service.
        """
        self.mdl = None

    def api(self, api_name, input):
        """
        Call an API with the given name and input.

        Parameters
        ----------
        api_name : str
            Name of the API to call.
        input : dict
            Input data for the API.

        Returns
        -------
        dict
            Response from the API.
        """
        method = getattr(self.mdl, api_name)
        return method(input)


class DummyService(BaseServing):
    """
    Dummy service class for testing purposes.

    Parameters
    ----------
    model_id : str
        The ID of the model to be served.
    config_json : dict, optional
        Configuration settings in JSON format.
    preferred_port : int, optional
        Preferred port for serving the model.
    url : str, optional
        URL for the served model.
    """

    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        BaseServing.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        Returns
        -------
        DummyService
            The instance of the service.
        """
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Exit the runtime context related to this object.

        Parameters
        ----------
        exception_type : type
            The exception type.
        exception_value : Exception
            The exception instance.
        traceback : traceback
        """
        self.close()

    def is_available(self):
        """
        Check if the dummy service is available.

        Returns
        -------
        bool
            True if the service is available, False otherwise.
        """
        return True

    def serve(self):
        """
        Serve the model using the dummy service.
        """
        pass

    def close(self):
        """
        Close the dummy service.
        """
        pass

    def api(self, api_name, input):
        """
        Call an API with the given name and input.

        Parameters
        ----------
        api_name : str
            Name of the API to call.
        input : dict
            Input data for the API.

        Returns
        -------
        dict
            Response from the API.
        """
        return self._api_with_url(api_name, input)


class PulledDockerImageService(BaseServing):
    """
    Service class for managing Docker image models pulled from a registry.

    Parameters
    ----------
    model_id : str
        The ID of the model to be served.
    config_json : dict, optional
        Configuration settings in JSON format.
    preferred_port : int, optional
        Preferred port for serving the model.
    url : str, optional
        URL for the served model.
    """

    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        self._is_docker_active()
        BaseServing.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )
        set_docker_host()
        self.client = docker.from_env()
        if preferred_port is None:
            self.port = find_free_port()
        else:
            self.port = preferred_port
        bundle_path = self._model_path(model_id)
        self.docker_tag = model_image_version_reader(bundle_path)
        self.logger.debug("Using port {0}".format(self.port))
        self.image_name = "{0}/{1}:{2}".format(
            DOCKERHUB_ORG, self.model_id, self.docker_tag
        )
        self.logger.debug("Starting Docker Daemon service")

        self.simple_docker = SimpleDocker()
        self.pid = -1
        self._mem_gb = self._get_memory()

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        Returns
        -------
        PulledDockerImageService
            The instance of the service.
        """
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Exit the runtime context related to this object.

        Parameters
        ----------
        exception_type : type
            The exception type.
        exception_value : Exception
            The exception instance.
        traceback : traceback
        """
        self.close()

    @throw_ersilia_exception()
    def _is_docker_active(self):
        dr = DockerRequirement()
        if not dr.is_active():
            raise DockerNotActiveError()
        return True

    def _api_with_url(self, api_name, input):
        if self.url is None:
            return
        self.logger.debug("Using URL: {0}".format(self.url))
        response = requests.post("{0}/{1}".format(self.url, api_name), json=input)
        return response.json()

    def _get_memory(self):
        info_file = os.path.join(self._model_path(self.model_id), INFORMATION_FILE)
        if not os.path.exists(info_file):
            return None
        with open(info_file, "r") as f:
            info = json.load(f)["card"]
        memory_field = "Memory Gb"
        if memory_field not in info:
            return None
        mem_gb = info[memory_field]
        if type(mem_gb) != int:
            return None
        if mem_gb < 2:
            mem_gb = 2
        self.logger.debug("Asking for {0} GB of memory".format(mem_gb))
        return mem_gb

    def is_available(self) -> bool:
        """
        Check if the Docker image service is available.

        Returns
        -------
        bool
            True if the service is available, False otherwise.
        """
        is_available = self.simple_docker.exists(
            DOCKERHUB_ORG, self.model_id, self.docker_tag
        )
        if is_available:
            self.logger.debug("Image {0} is available locally".format(self.image_name))
            return True
        else:
            self.logger.debug(
                "Image {0} is not available locally".format(self.image_name)
            )
            return False

    def _stop_all_containers_of_image(self):
        self.logger.debug(
            "Stopping all containers related to model {0}".format(self.model_id)
        )
        containers = self.client.containers.list(all=True)
        for container in containers:
            self.logger.debug(f"Checking container {container.name}")
            if container.name.startswith(self.model_id):
                self.logger.debug(
                    "Stopping and removing container {0}".format(container.name)
                )
                self._delete_temp_files(container.name)
                container.stop()
                self.logger.debug("Container stopped")
                container.remove()
                self.logger.debug("Container removed")

    @throw_ersilia_exception()
    def _get_apis(self):
        file_name = os.path.join(
            self._get_bundle_location(self.model_id), APIS_LIST_FILE
        )
        self.logger.debug("Getting APIs")
        if os.path.exists(file_name):
            with open(file_name, "r") as f:
                apis_list = []
                for l in f:
                    apis_list += [l.rstrip()]
            if len(apis_list) > 0:
                return apis_list

        def _url_exists(url):
            try:
                response = requests.head(url, allow_redirects=True)
                if response.status_code in [200, 301, 302]:
                    self.logger.debug(f"The URL {url} exists.")
                    return True
                else:
                    self.logger.debug(
                        f"The URL {url} does not exist. Status code: {response.status_code}"
                    )
                    return False

            except requests.exceptions.RequestException as e:
                self.logger.debug(f"An error occurred: {e}")
                return False

        self.logger.debug("Trying to get them using info endpoint")
        if _url_exists(f"{self.url}/info"):
            url = "{0}/info".format(self.url)
            self.logger.debug("Using URL: {0}".format(url))
            data = "{}"
            response = requests.post(url, data=data)
            self.logger.debug("Status code: {0}".format(response.status_code))
            if response.status_code == 502:
                raise BadGatewayError(url)
            elif response.status_code == 405:
                response = requests.get(url)
            else:
                response.raise_for_status()
            apis_list = json.loads(response.text)["apis_list"]
        else:
            apis_list = []
            for api in ALLOWED_API_NAMES:
                github_base_url = "https://raw.githubusercontent.com/ersilia-os/{0}/refs/heads/main/model/framework".format(
                    self.model_id
                )
                github_url = "{0}/{1}.sh".format(github_base_url, api)
                self.logger.debug("Checking URL: {0}".format(github_url))
                response = requests.head(github_url)
                if response.status_code == 200:
                    apis_list.append(api)

        self.logger.debug("Writing file {0}".format(file_name))
        with open(file_name, "w") as f:
            for api in apis_list:
                f.write(api + os.linesep)
        return apis_list

    def is_url_available(self, url):
        """
        Check if the given URL is available.

        Parameters
        ----------
        url : str
            The URL to check.

        Returns
        -------
        bool
            True if the URL is available, False otherwise.
        """
        try:
            response = requests.get(url, timeout=5)
            response.raise_for_status()
        except requests.HTTPError:
            return False
        except Exception:
            return False
        else:
            return True

    def _wait_until_container_is_running(self):
        while True:
            if self.is_url_available(self.url):
                break
            else:
                self.logger.debug("Container in {0} is not ready yet".format(self.url))
            time.sleep(1)

    def _create_docker_network(self):
        existing_networks = self.client.networks.list(
            filters={"name": DEFAULT_DOCKER_NETWORK_NAME}
        )
        if existing_networks:
            self.network = existing_networks[0]
            self.logger.info(f"Docker network already exists: {self.network.name}")
        else:
            self.network = self.client.networks.create(
                name=DEFAULT_DOCKER_NETWORK_NAME,
                driver=DEFAULT_DOCKER_NETWORK_BRIDGE,
                attachable=True,
            )
            self.logger.info(f"Docker network has been created: {self.network.name}")

    def serve(self):
        """
        Serve the model using the Docker image service.
        """
        self._create_docker_network()
        self._stop_all_containers_of_image()
        self.container_name = f"{self.model_id}_{str(uuid.uuid4())[:4]}"

        env = {
            "REDIS_HOST": os.getenv("REDIS_HOST", "redis"),
            "REDIS_PORT": os.getenv("REDIS_PORT", "6379"),
            "REDIS_URI": os.getenv("REDIS_URI", "redis://redis:6379"),
            "REDIS_EXPIRATION": os.getenv("REDIS_EXPIRATION", str(3600 * 24 * 7)),
        }

        run_kwargs = dict(
            image=self.image_name,
            name=self.container_name,
            detach=True,
            ports={"80/tcp": self.port},
            environment=env,
            network=DEFAULT_DOCKER_NETWORK_NAME,
        )

        if self._mem_gb is not None:
            run_kwargs["mem_limit"] = f"{self._mem_gb}g"

        self.logger.debug(f"Running container with env: {env!r}")
        self.container = self.client.containers.run(**run_kwargs)

        self.container_id = self.container.id
        self.url = f"http://0.0.0.0:{self.port}"
        self._wait_until_container_is_running()
        self._apis_list = self._get_apis()
        self.logger.debug(self._apis_list)

    def api(self, api_name: str, input: dict) -> dict:
        """
        Call an API with the given name and input.

        Parameters
        ----------
        api_name : str
            Name of the API to call.
        input : dict
            Input data for the API.

        Returns
        -------
        dict
            Response from the API.
        """
        return self._api_with_url(api_name=api_name, input=input)

    def _delete_temp_files(self, container_id):
        """
        Deletes a file inside a running Docker container.

        :param container_id: str - The ID or name of the container.
        :param file_path: str - The absolute path of the file to delete inside the container.
        :return: None
        """
        # Create a Docker client from environment variables
        set_docker_host()
        client = docker.from_env()
        container = client.containers.get(container_id)

        # Create the command to remove the file.
        # Using 'rm -rf' avoids errors if the file doesn't exist.
        ls_cmd = "ls /tmp"
        rem_cmd = "rm -rf /tmp/{0}"

        try:
            # Execute the command inside the container
            exec_result = container.exec_run(ls_cmd)
            # Check if the command executed successfully
            if exec_result.exit_code == 0:
                list_dirs = exec_result.output.decode("utf-8").split("\n")
                for f in list_dirs:
                    container.exec_run(rem_cmd.format(f))
                    self.logger.debug(
                        f"Deleted temp file {f} from container {container_id}"
                    )
            else:
                output = (
                    exec_result.output.decode("utf-8")
                    if exec_result.output
                    else "No output"
                )
                self.logger.debug(
                    f"Error deleting file. Exit code: {exec_result.exit_code}. Output: {output}"
                )
        except Exception as e:
            self.logger.debug(
                f"An error occurred during execution inside the container: {e}"
            )

    def close(self):
        """
        Close the Docker image service by stopping and removing the container.
        """
        self.logger.debug("Stopping and removing container")
        # Here we remove temp files so that they are not left behind
        self._stop_all_containers_of_image()


class HostedService(BaseServing):
    """
    Service class for managing hosted models.

    Parameters
    ----------
    model_id : str
        The ID of the model to be served.
    config_json : dict, optional
        Configuration settings in JSON format.
    preferred_port : int, optional
        Preferred port for serving the model.
    url : str, optional
        URL for the served model.
    """

    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        BaseServing.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=None,
        )
        if url is None:
            self.url = self._resolve_url()
        else:
            self.url = url
        self.pid = -1

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        Returns
        -------
        HostedService
            The instance of the service.
        """
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Exit the runtime context related to this object.

        Parameters
        ----------
        exception_type : type
            The exception type.
        exception_value : Exception
            The exception instance.
        traceback : traceback
        """
        self.close()

    def _api_with_url(self, api_name, input):
        if self.url is None:
            return
        self.logger.debug("Using URL: {0}".format(self.url))
        response = requests.post("{0}/{1}".format(self.url, api_name), json=input)
        return response.json()

    def _resolve_url(self):
        from_hosted_file = os.path.join(
            self._model_path(self.model_id), IS_FETCHED_FROM_HOSTED_FILE
        )
        self.logger.debug("Reading hosted file: {0}".format(from_hosted_file))
        if not os.path.exists(from_hosted_file):
            return None
        with open(from_hosted_file, "r") as f:
            data = json.load(f)
        self.logger.debug("From hosted file: {0}".format(data))
        if not data["hosted"]:
            return None
        return data["url"]

    def is_available(self) -> bool:
        """
        Check if the hosted service is available.

        Returns
        -------
        bool
            True if the service is available, False otherwise.
        """
        if self.is_url_available(self.url):
            self.logger.debug("URL {0} is available".format(self.url))
            return True
        else:
            self.logger.debug("URL {0} is not available".format(self.url))
            return False

    def _get_apis(self):
        """
        Get the list of APIs available for the model.

        Returns
        -------
        list
            List of API names.
        """
        file_name = os.path.join(
            self._get_bundle_location(self.model_id), APIS_LIST_FILE
        )
        self.logger.debug("Getting APIs")
        if os.path.exists(file_name):
            with open(file_name, "r") as f:
                apis_list = []
                for l in f:
                    apis_list += [l.rstrip()]
            if len(apis_list) > 0:
                return apis_list

        def _url_exists(url):
            try:
                response = requests.head(url, allow_redirects=True)
                if response.status_code in [200, 301, 302]:
                    self.logger.debug(f"The URL {url} exists.")
                    return True
                else:
                    self.logger.debug(
                        f"The URL {url} does not exist. Status code: {response.status_code}"
                    )
                    return False

            except requests.exceptions.RequestException as e:
                self.logger.debug(f"An error occurred: {e}")
                return False

        self.logger.debug("Trying to get them using info endpoint")
        if _url_exists(f"{self.url}/info"):
            url = "{0}/info".format(self.url)
            self.logger.debug("Using URL: {0}".format(url))
            data = "{}"
            response = requests.post(url, data=data)
            self.logger.debug("Status code: {0}".format(response.status_code))
            if response.status_code == 502:
                raise BadGatewayError(url)
            elif response.status_code == 405:
                response = requests.get(url)
            else:
                response.raise_for_status()
            apis_list = json.loads(response.text)["apis_list"]
        else:
            apis_list = []
            for api in ALLOWED_API_NAMES:
                github_base_url = "https://raw.githubusercontent.com/ersilia-os/{0}/refs/heads/main/model/framework".format(
                    self.model_id
                )
                github_url = "{0}/{1}.sh".format(github_base_url, api)
                self.logger.debug("Checking URL: {0}".format(github_url))
                response = requests.head(github_url)
                if response.status_code == 200:
                    apis_list.append(api)

        self.logger.debug("Writing file {0}".format(file_name))
        with open(file_name, "w") as f:
            for api in apis_list:
                f.write(api + os.linesep)
        return apis_list

    def is_url_available(self, url):
        """
        Check if the given URL is available.

        Parameters
        ----------
        url : str
            The URL to check.

        Returns
        -------
        bool
            True if the URL is available, False otherwise.
        """
        try:
            response = requests.get(url, timeout=5)
            response.raise_for_status()
        except requests.HTTPError:
            return False
        except Exception:
            return False
        else:
            return True

    def serve(self):
        """
        Serve the model using the hosted service.
        """
        self._apis_list = self._get_apis()
        self.logger.debug(self._apis_list)

    def api(self, api_name: str, input: dict) -> dict:
        """
        Call an API with the given name and input.

        Parameters
        ----------
        api_name : str
            Name of the API to call.
        input : dict
            Input data for the API.

        Returns
        -------
        dict
            Response from the API.
        """
        return self._api_with_url(api_name=api_name, input=input)

    def close(self):
        """
        Close the hosted service.
        """
        pass
