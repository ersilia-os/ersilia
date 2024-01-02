import tempfile
import os
import json
import time
import importlib
import requests
import uuid
import docker
from .. import ErsiliaBase, throw_ersilia_exception
from ..utils.terminal import run_command
from ..utils.ports import find_free_port
from ..db.environments.localdb import EnvironmentDb
from ..db.environments.managers import DockerManager
from ..utils.conda import SimpleConda
from ..utils.docker import SimpleDocker
from ..utils.venv import SimpleVenv
from ..default import DEFAULT_VENV
from ..default import PACKMODE_FILE, APIS_LIST_FILE
from ..default import DOCKERHUB_ORG, DOCKERHUB_LATEST_TAG
from ..default import IS_FETCHED_FROM_HOSTED_FILE
from ..default import INFORMATION_FILE
from ..utils.exceptions_utils.serve_exceptions import BadGatewayError

SLEEP_SECONDS = 1
TIMEOUT_SECONDS = 1000


class BaseServing(ErsiliaBase):
    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.bundle_tag = self._get_latest_bundle_tag(model_id=self.model_id)
        self.port = preferred_port

    def _get_info_from_bento(self):
        """Get info available from the Bento"""
        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "info.json")
        cmd = "bentoml info --quiet {0}:{1} > {2}".format(
            self.model_id, self.bundle_tag, tmp_file
        )
        self.logger.debug(
            "Getting info from BentoML and storing in {0}".format(tmp_file)
        )
        run_command(cmd)
        with open(tmp_file, "r") as f:
            info = json.load(f)
        self.logger.debug("Info {0}".format(info))
        return info

    def _get_apis_from_bento(self):
        """Get APIs available for the model, according to the info Bento"""
        self.logger.debug("Getting APIs from Bento")
        info = self._get_info_from_bento()
        apis_list = []
        for item in info["apis"]:
            apis_list += [item["name"]]
        return apis_list

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
        print(apis_list)
        if len(apis_list) > 0:
            return apis_list
        else:
            return None

    def _get_apis_from_where_available(self):
        apis_list = self._get_apis_from_apis_list()
        if apis_list is None:
            apis_list = self._get_apis_from_bento()
        if apis_list is None:
            apis_list = []
        for api in apis_list:
            yield api

    def _api_with_url(self, api_name, input):
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

    def _bentoml_serve(self, runcommand_func=None):
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
        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
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

    def _close(self):
        try:
            os.kill(self.pid, 9)
        except:
            self.logger.info("PID {0} is unassigned".format(self.pid))


class SystemBundleService(_BentoMLService):
    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        _BentoMLService.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def _run_command(self, cmd):
        return run_command(cmd)

    def is_available(self):
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
        self._bentoml_serve()

    def close(self):
        self._close()

    def api(self, api_name, input):
        return self._api_with_url(api_name, input)


class VenvEnvironmentService(_BentoMLService):
    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        _BentoMLService.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )
        self.venv = SimpleVenv(self._model_path(model_id))

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def _model_path(self, model_id):
        return os.path.join(self._dest_dir, model_id)

    def _run_command(self, cmd):
        return self.venv.run_commandlines(DEFAULT_VENV, cmd)

    def is_available(self):
        if not self.venv.exists(DEFAULT_VENV):
            return False

    def serve(self):
        self._bentoml_serve(self._run_command)

    def close(self):
        self._close()

    def api(self, api_name, input):
        return self._api_with_url(api_name, input)


class CondaEnvironmentService(_BentoMLService):
    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        _BentoMLService.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )
        self.db = EnvironmentDb()
        self.db.table = "conda"
        self.conda = SimpleConda()

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def _get_env_name(self):
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
        env = self._get_env_name()
        if env is not None:
            return True
        else:
            return False

    def serve(self):
        self._bentoml_serve(self._run_command)

    def close(self):
        self._close()

    def api(self, api_name, input):
        return self._api_with_url(api_name, input)


class DockerImageService(BaseServing):
    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
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
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def _get_env_name(self):
        envs = list(self.db.envs_of_model(self.model_id))
        for env in envs:
            org, img, tag = self.docker._splitter(env)
            if self.docker.exists(org=org, img=img, tag=tag):
                self.logger.debug("Docker image found {0}".format(env))
                return env
        return None

    def is_available(self):
        env = self._get_env_name()
        if env is not None:
            self.logger.debug("Docker image service available")
            return True
        else:
            self.logger.debug("Docker image service not available")
            return False

    def serve(self):
        self.logger.debug("Calling docker manager")
        res = self.dm.run(self.model_id)
        self.container_name = res["container_name"]
        self.port = res["port"]
        self.url = "http://0.0.0.0:{0}".format(self.port)

    def close(self):
        self.df.stop_containers(self.model_id)

    def api(self, api_name, input):
        return self._api_with_url(api_name, input)


# TODO: Include 'pip' within available service_class
class PipInstalledService(BaseServing):
    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        BaseServing.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )
        self.pid = -1

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def _import(self):
        try:
            model = importlib.import_module(self.model_id, package=None)
            return model
        except ModuleNotFoundError:
            return None

    def is_available(self):
        model = self._import()
        if model is not None:
            return True
        else:
            return False

    def serve(self):
        model = self._import()
        self.mdl = model.load()

    def close(self):
        self.mdl = None

    def api(self, api_name, input):
        method = getattr(self.mdl, api_name)
        return method(input)


class DummyService(BaseServing):
    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        BaseServing.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def is_available(self):
        return True

    def serve(self):
        pass

    def close(self):
        pass

    def api(self, api_name, input):
        return self._api_with_url(api_name, input)


class PulledDockerImageService(BaseServing):
    def __init__(self, model_id, config_json=None, preferred_port=None, url=None):
        BaseServing.__init__(
            self,
            model_id=model_id,
            config_json=config_json,
            preferred_port=preferred_port,
        )
        self.client = docker.from_env()
        if preferred_port is None:
            self.port = find_free_port()
        else:
            self.port = preferred_port
        self.logger.debug("Using port {0}".format(self.port))
        self.image_name = "{0}/{1}:{2}".format(
            DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG
        )
        self.logger.debug("Starting Docker Daemon service")
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        self.logger.debug(
            "Creating temporary folder {0} and mounting as volume in container".format(
                self.tmp_folder
            )
        )
        self.simple_docker = SimpleDocker()
        self.pid = -1
        self._mem_gb = self._get_memory()

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

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

    def is_available(self):
        is_available = self.simple_docker.exists(
            DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG
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
            if container.name.startswith(self.model_id):
                self.logger.debug(
                    "Stopping and removing container {0}".format(container.name)
                )
                container.stop()
                self.logger.debug("Container stopped")
                container.remove()
                self.logger.debug("Container removed")

    @throw_ersilia_exception
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
        self.logger.debug("Getting them using info endpoint")
        url = "{0}/info".format(self.url)
        self.logger.debug("Using URL: {0}".format(url))
        data = "{}"
        response = requests.post(url, data=data)
        self.logger.debug("Status code: {0}".format(response.status_code))
        if response.status_code == 502:
            raise BadGatewayError(url)
        apis_list = json.loads(response.text)["apis_list"]
        self.logger.debug("Writing file {0}".format(file_name))
        with open(file_name, "w") as f:
            for api in apis_list:
                f.write(api + os.linesep)
        return apis_list

    def is_url_available(self, url):
        try:
            response = requests.get(url, timeout=5)
            response.raise_for_status()
        except requests.HTTPError as http_err:
            return False
        except Exception as err:
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

    def serve(self):
        self._stop_all_containers_of_image()
        self.container_name = "{0}_{1}".format(self.model_id, str(uuid.uuid4())[:4])
        self.volumes = {self.tmp_folder: {"bind": "/ersilia_tmp", "mode": "rw"}}
        self.logger.debug("Trying to run container")
        if self._mem_gb is None:
            self.container = self.client.containers.run(
                self.image_name,
                name=self.container_name,
                detach=True,
                ports={"80/tcp": self.port},
                volumes=self.volumes,
            )
        else:
            self.container = self.client.containers.run(
                self.image_name,
                name=self.container_name,
                detach=True,
                ports={"80/tcp": self.port},
                volumes=self.volumes,
                mem_limit="{0}g".format(self._mem_gb),
            )
        self.logger.debug("Serving container {0}".format(self.container_name))
        self.container_id = self.container.id
        self.logger.debug("Running container {0}".format(self.container_id))
        self.url = "http://0.0.0.0:{0}".format(self.port)
        self._wait_until_container_is_running()
        self._apis_list = self._get_apis()
        self.logger.debug(self._apis_list)

    def api(self, api_name, input):
        return self._api_with_url(api_name=api_name, input=input)

    def close(self):
        self.logger.debug("Stopping and removing container")
        self._stop_all_containers_of_image()


class HostedService(BaseServing):
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
        return self

    def __exit__(self, exception_type, exception_value, traceback):
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

    def is_available(self):
        if self.is_url_available(self.url):
            self.logger.debug("URL {0} is available".format(self.url))
            return True
        else:
            self.logger.debug("URL {0} is not available".format(self.url))
            return False

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
        self.logger.debug("Getting them using info endpoint")
        url = "{0}/info".format(self.url)
        self.logger.debug("Using URL: {0}".format(url))
        data = "{}"
        apis_list = json.loads(requests.post(url, data=data).text)["apis_list"]
        self.logger.debug("Writing file {0}".format(file_name))
        with open(file_name, "w") as f:
            for api in apis_list:
                f.write(api + os.linesep)
        return apis_list

    def is_url_available(self, url):
        try:
            response = requests.get(url, timeout=5)
            response.raise_for_status()
        except requests.HTTPError as http_err:
            return False
        except Exception as err:
            return False
        else:
            return True

    def serve(self):
        self._apis_list = self._get_apis()
        self.logger.debug(self._apis_list)

    def api(self, api_name, input):
        return self._api_with_url(api_name=api_name, input=input)

    def close(self):
        pass
