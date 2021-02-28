import tempfile
import os
import json
import time
import subprocess
import importlib
from .. import ErsiliaBase
from ..utils.terminal import run_command
from ..utils.ports import find_free_port

SLEEP_SECONDS = 1


class SystemBundleService(ErsiliaBase):

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.bundle_tag = self._get_latest_bundle_tag(model_id=self.model_id)

    def _serve(self):
        """Simply try to serve model with bentoml, locally"""
        tmp_folder = tempfile.mkdtemp()
        tmp_folder = "/Users/mduran/Desktop/mytest" # remove
        if not os.path.exists(tmp_folder): os.mkdir(tmp_folder) # remove
        tmp_script = os.path.join(tmp_folder, "serve.sh")
        tmp_file = os.path.join(tmp_folder, "serve.log")
        tmp_pid = os.path.join(tmp_folder, "serve.pid")
        sl = ['#!/bin/bash']
        sl += ['bentoml serve {0}:{1} --port {2} &> {3} &'.format(self.model_id, self.bundle_tag, self.port(), tmp_file)]
        sl += ['_pid=$!']
        sl += ['echo "$_pid" > {0}'.format(tmp_pid)]
        with open(tmp_script, "w") as f:
            for l in sl:
                f.write(l+os.linesep)
        cmd = "bash {0}".format(tmp_script)
        run_command(cmd, quiet=True)
        with open(tmp_pid, "r") as f:
            pid = int(f.read().strip())
        while not os.path.exists(tmp_file):
            time.sleep(SLEEP_SECONDS)
        with open(tmp_file, "r") as f:
            print(f.read())

    def is_available(self):
        # TODO
        pass

    def serve(self):
        # TODO
        pass

    def close(self):
        # TODO
        pass


class CondaEnvironmentService(ErsiliaBase):

    from ..db.environments.localdb import EnvironmentDb
    from ..utils.conda import SimpleConda

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.db = EnvironmentDb()
        self.conda = SimpleConda()

    def _get_env_name(self):
        envs = list(self.db.envs_of_model(self.model_id))
        for env in envs:
            if self.conda.exists(env):
                return env
        return None

    def is_available(self):
        env = self._get_env_name()
        if env is not None:
            return True
        else:
            return False

    def serve(self):
        # TODO
        pass

    def close(self):
        # TODO
        pass


class DockerImageService(ErsiliaBase):

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id

    def is_available(self):
        pass

    def serve(self):
        pass

    def close(self):
        pass


class PipInstalledService(ErsiliaBase):

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id

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


class ServingModality(ErsiliaBase):

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, model_id=model_id, config_json=config_json)

    def _is_remote(self):
        pass

    def _is_bento(self):
        pass

    def _is_pip(self):
        pass

    def _is_conda(self):
        pass

    def _is_docker(self):
        pass

    def _get_info(self):
        """Get info available from the Bento"""
        tmp_folder = tempfile.mkdtemp()
        tmp_file = os.path.join(tmp_folder, "info.json")
        cmd = "bentoml info --quiet {0}:{1} > {2}".format(self.model_id, self.bundle_tag, tmp_file)
        run_command(cmd, quiet=True)
        with open(tmp_file, "r") as f:
            info = json.load(f)
        return info

    def _get_apis(self):
        """Get APIs available for the model, according to the info Bento"""
        info = self._get_info()
        for item in info["apis"]:
            yield item["name"]
