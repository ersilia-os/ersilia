import tempfile
import os
import json
import time
import subprocess
import importlib
import requests
from .. import ErsiliaBase
from ..utils.terminal import run_command
from ..utils.ports import find_free_port
from ..db.environments.localdb import EnvironmentDb
from ..db.environments.managers import DockerManager
from ..utils.conda import SimpleConda
from ..utils.docker import SimpleDocker
from ..default import DOCKERHUB_ORG

SLEEP_SECONDS = 1
TIMEOUT_SECONDS = 1000


class BaseServing(ErsiliaBase):

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.bundle_tag = self._get_latest_bundle_tag(model_id=self.model_id)

    def _get_info_from_bento(self):
        """Get info available from the Bento"""
        tmp_folder = tempfile.mkdtemp()
        tmp_file = os.path.join(tmp_folder, "info.json")
        cmd = "bentoml info --quiet {0}:{1} > {2}".format(self.model_id, self.bundle_tag, tmp_file)
        run_command(cmd, quiet=True)
        with open(tmp_file, "r") as f:
            info = json.load(f)
        return info

    def _get_apis_from_bento(self):
        """Get APIs available for the model, according to the info Bento"""
        info = self._get_info_from_bento()
        for item in info["apis"]:
            yield item["name"]

    def _api_with_url(self, api_name, input):
        if self.url is None:
            return
        response = requests.post("{0}/{1}".format(self.url, api_name), json=input)
        return response.json()


class _BentoMLService(BaseServing):

    def __init__(self, model_id, config_json=None):
        BaseServing.__init__(self, model_id=model_id, config_json=config_json)
        self.SEARCH_PRE_STRING = "* Running on "
        self.SEARCH_SUF_STRING = "(Press CTRL+C to quit)"
        self.ERROR_STRING = "error"

    def _bentoml_serve(self, runcommand_func):
        """Simply try to serve model with bentoml, locally"""
        self.port = find_free_port()
        tmp_folder = tempfile.mkdtemp()
        tmp_script = os.path.join(tmp_folder, "serve.sh")
        tmp_file = os.path.join(tmp_folder, "serve.log")
        tmp_pid = os.path.join(tmp_folder, "serve.pid")
        sl = ['#!/bin/bash']
        sl += ['bentoml serve {0}:{1} --port {2} &> {3} &'.format(self.model_id, self.bundle_tag, self.port, tmp_file)]
        sl += ['_pid=$!']
        sl += ['echo "$_pid" > {0}'.format(tmp_pid)]
        with open(tmp_script, "w") as f:
            for l in sl:
                f.write(l+os.linesep)
        cmd = "bash {0}".format(tmp_script)
        run_command(cmd, quiet=True)
        with open(tmp_pid, "r") as f:
            self.pid = int(f.read().strip())
        for _ in range(int(TIMEOUT_SECONDS/SLEEP_SECONDS)):
            # If error string is identified, finish
            with open(tmp_file, "r") as f:
                r = f.read()
                if self.ERROR_STRING in r:
                    self.url = None
                    return
            # If everything looks good, wait until server is ready
            with open(tmp_file, "r") as f:
                r = f.read()
                if self.SEARCH_PRE_STRING not in r or self.SEARCH_SUF_STRING not in r:
                    time.sleep(SLEEP_SECONDS)
                    continue
            # When the search strings are found get url
            with open(tmp_file, "r") as f:
                for l in f:
                    if self.SEARCH_PRE_STRING in l and self.SEARCH_SUF_STRING in l:
                        self.url = l.split(self.SEARCH_PRE_STRING)[1].split(" ")[0]
                        return
        self.url = None

    def _close(self):
        cmd = "kill {0}".format(self.pid)
        run_command(cmd, quiet=True)


class SystemBundleService(_BentoMLService):

    def __init__(self, model_id, config_json=None):
        _BentoMLService.__init__(self, model_id=model_id, config_json=config_json)

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def _serve(self):
        self._bentoml_serve(run_command)

    def is_available(self):
        self._serve()
        self.close()
        if self.url is not None:
            avail = True
        else:
            avail = False
        self.url = None
        return avail

    def serve(self):
        self._serve()

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

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    @staticmethod
    def _splitter(env):
        return env.split(":")

    def _get_env_name(self):
        envs = list(self.db.envs_of_model(self.model_id))
        for env in envs:
            img, tag = self._splitter(env)
            if self.docker.exists(org=DOCKERHUB_ORG, img=img, tag=tag):
                return env
        return None

    def is_available(self):
        env = self._get_env_name()
        if env is not None:
            return True
        else:
            return False

    def serve(self):
        res = self.dm.run(self.model_id)
        self.container_name = res["container_name"]
        self.port = res["port"]
        self.url = "http://localhost:{0}".format(self.port)

    def close(self):
        self.dm._delete_container(self.container_name)

    def api(self, api_name, input):
        return self._api_with_url(api_name, input)


class PipInstalledService(BaseServing):

    def __init__(self, model_id, config_json=None):
        BaseServing.__init__(self, model_id=model_id, config_json=config_json)

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
