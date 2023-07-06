import os
import importlib
import json

from ...utils.docker import SimpleDocker
from ...utils.conda import SimpleConda
from ...db.environments.localdb import EnvironmentDb
from ... import ErsiliaBase
from ...default import IS_FETCHED_FROM_DOCKERHUB_FILE


class ModelStatus(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def is_downloaded(self, model_id):
        essentials = ["README.md", "model"]  #  essential files
        dst_dir = os.path.join(self._dest_dir, model_id)
        if not os.path.exists(dst_dir):
            return False
        items = {i for i in os.listdir(dst_dir)}
        for essential in essentials:
            if essential not in items:
                return False
        return True

    def is_docker(self, model_id):
        db = EnvironmentDb(config_json=self.config_json)
        db.table = "docker"
        docker = SimpleDocker()
        envs = db.envs_of_model(model_id)
        for env in envs:
            img, tag = env.split(":")
            if docker.exists(self.cfg.EXT.DOCKERHUB_ORG, img, tag):
                return True
        return False

    def is_pulled_docker(self, model_id):
        model_dir = os.path.join(self._model_path(model_id=model_id))
        json_file = os.path.join(model_dir, IS_FETCHED_FROM_DOCKERHUB_FILE)
        if not os.path.exists(json_file):
            return False
        with open(json_file, "r") as f:
            data = json.load(f)
        return data["docker_hub"]

    def is_conda(self, model_id):
        db = EnvironmentDb(config_json=self.config_json)
        db.table = "conda"
        conda = SimpleConda()
        envs = db.envs_of_model(model_id)
        for env in envs:
            if conda.exists(env):
                return True
        return False

    def is_pip(self, model_id):
        try:
            importlib.import_module(model_id, package=None)
            return True
        except ModuleNotFoundError:
            return False

    def _is_bento_folder(self, model_folder):
        if model_folder is None:
            return False
        if not os.path.exists(model_folder):
            return False
        essentials = ["bentoml.yml"]
        items = {i for i in os.listdir(model_folder)}
        for essential in essentials:
            if essential not in items:
                return False
        return True

    def is_bundle(self, model_id):
        model_folder = self._get_bundle_location(model_id)
        return self._is_bento_folder(model_folder)

    def is_bentoml(self, model_id):
        model_folder = self._get_bentoml_location(model_id)
        return self._is_bento_folder(model_folder)

    def status(self, model_id):
        """Check installation and deployment status of the model"""
        results = {
            "download": self.is_downloaded(model_id),
            "bentoml": self.is_bentoml(model_id),
            "bundle": self.is_bundle(model_id),
            "docker": self.is_docker(model_id),
            "pulled_docker": self.is_pulled_docker(model_id),
            "conda": self.is_conda(model_id),
            "pip": self.is_pip(model_id),
        }
        return results
