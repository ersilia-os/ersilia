import tempfile
import os
import json
import time
import importlib
import requests
from .. import ErsiliaBase
from ..utils.terminal import run_command
from ..utils.ports import find_free_port
from ..db.environments.localdb import EnvironmentDb
from ..db.environments.managers import DockerManager
from ..utils.conda import SimpleConda
from ..utils.docker import SimpleDocker
from ..utils.venv import SimpleVenv
from ..default import DEFAULT_VENV
from ..default import PACKMODE_FILE

SLEEP_SECONDS = 1
TIMEOUT_SECONDS = 1000


class BaseServing(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.bundle_tag = self._get_latest_bundle_tag(model_id=self.model_id)

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
        for item in info["apis"]:
            yield item["name"]

    def _api_with_url(self, api_name, input):
        if self.url is None:
            return
        self.logger.debug("Using URL: {0}".format(self.url))
        response = requests.post("{0}/{1}".format(self.url, api_name), json=input)
        return response.json()


class _BentoMLService(BaseServing):
    def __init__(self, model_id, config_json=None):
        BaseServing.__init__(self, model_id=model_id, config_json=config_json)
        self.SEARCH_PRE_STRING = "* Running on "
        self.SEARCH_SUF_STRING = "Press CTRL+C to quit"
        self.ERROR_STRING = "error"

    def _bentoml_serve(self, runcommand_func=None):
        self.logger.debug("Trying to serve model with BentoML locally")
        self.port = find_free_port()
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
    def __init__(self, model_id, config_json=None):
        _BentoMLService.__init__(self, model_id=model_id, config_json=config_json)

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
    def __init__(self, model_id, config_json=None):
        _BentoMLService.__init__(self, model_id=model_id, config_json=config_json)
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
    def __init__(self, model_id, config_json=None):
        _BentoMLService.__init__(self, model_id=model_id, config_json=config_json)
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
    def __init__(self, model_id, config_json=None):
        BaseServing.__init__(self, model_id=model_id, config_json=config_json)
        self.db = EnvironmentDb()
        self.db.table = "docker"
        self.docker = SimpleDocker()
        self.dm = DockerManager(config_json=config_json)
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
    def __init__(self, model_id, config_json=None):
        BaseServing.__init__(self, model_id=model_id, config_json=config_json)
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
    def __init__(self, model_id, config_json=None):
        BaseServing.__init__(self, model_id=model_id, config_json=config_json)

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
