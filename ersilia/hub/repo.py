import os
import json
from .. import ErsiliaBase
from ..utils.paths import model_id_from_path


class RepoUtils(ErsiliaBase):

    def __init__(self, path, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.dockerhub_org = self.cfg.EXT.DOCKERHUB_ORG
        self.config_in_img = os.path.join(self.cfg.HUB.IMAGE_REPODIR, self.cfg.HUB.CONFIG_FILE)
        self.path = os.path.normpath(os.path.dirname(os.path.abspath(path)))

    def _root_path(self):
        model_id = self._get_model_id_from_path()
        if model_id is None:
            return self.path
        else:
            splits = self.path.split(os.sep)
            idx = splits.index(model_id)
            idx += 1
            return os.sep.join(splits[:idx])

    def _get_model_id_from_path(self):
        return model_id_from_path(self.path)

    def _get_model_id_from_config(self):
        if not os.path.exists(self.config_in_img):
            return None
        with open(self.config_in_img, "r") as f:
            cfg = json.load(f)
        return cfg["model_id"].replace("'", "").replace('"', '')

    def get_model_id(self):
        model_id = self._get_model_id_from_path()
        if model_id is None:
            return self._get_model_id_from_config()
        else:
            return model_id

    def get_conda_env_yml_file(self):
        root = self._root_path()
        return os.path.join(root, "environment.yml")

    def get_docker_repo_image(self, model_id):
        return os.path.join(self.dockerhub_org, "%s:repo" % model_id)

    @staticmethod
    def rename_service(model_id):
        cmd = "Service.__name__ = '%s'\n" % model_id
        cmd += "%s = Service" % model_id
        return cmd
