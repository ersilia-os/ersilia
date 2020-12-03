import os
import json
from .. import ErsiliaBase
from ..utils.paths import model_id_from_path
from ..default import CONDA_ENV_YML_FILE

ROOT_CHECKFILE = "README.md"


class RepoUtils(ErsiliaBase):

    def __init__(self, path, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.dockerhub_org = self.cfg.EXT.DOCKERHUB_ORG
        self.config_in_img = os.path.join(self.cfg.ENV.DOCKER.IMAGE_REPODIR, self.cfg.HUB.CONFIG_FILE)
        self.path = os.path.normpath(os.path.dirname(os.path.abspath(path)))

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

    def _root_path(self):
        model_id = self._get_model_id_from_path()
        if model_id is None:
            return self.path
        splits = self.path.split(os.sep)
        idx = splits.index(model_id)
        idx += 1
        root = os.sep.join(splits[:idx])
        # check the we really are in the root, and not in a version (as done by BentoML)
        files = os.listdir(root)
        if ROOT_CHECKFILE in files:
            return root
        # try to find bundles dir
        elif self._bundles_dir in root:
            tag = self._get_latest_bundle_tag(model_id)
            return os.path.join(root, tag)
        # or must be in the bentoml path
        elif self._bentoml_dir in root:
            tag = self._get_latest_bentoml_tag(model_id)
            return os.path.join(root, tag)
        else:
            return None

    def get_conda_env_yml_file(self):
        root = self._root_path()
        if root is None:
            model_id = self.get_model_id()
            # try to find yml in bundles
            yml = os.path.join(self._bundles_dir, model_id, self._get_latest_bundle_tag(model_id), CONDA_ENV_YML_FILE)
            if os.path.exists(yml):
                return yml
            # try to find yml in bentoml
            yml = os.path.join(self._bentoml_dir, model_id, self._get_latest_bentoml_tag(model_id), CONDA_ENV_YML_FILE)
            if os.path.exists(yml):
                return yml
        else:
            return os.path.join(root, CONDA_ENV_YML_FILE)


    def get_docker_repo_image(self, model_id):
        return os.path.join(self.dockerhub_org, "{0}:{1}".format(model_id, self.cfg.ENV.DOCKER.REPO_TAG))

    @staticmethod
    def rename_service(model_id):
        cmd = "Service.__name__ = '%s'\n" % model_id
        cmd += "%s = Service" % model_id
        return cmd
