import os
import shutil

from . import BaseAction
from ....utils.download import GitHubDownloader
from ....utils.paths import Paths
from ...bundle.repo import PackFile


MODEL_DIR = "model"


class ModelRepositoryGetter(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.token = self.cfg.HUB.TOKEN
        self.github_down = GitHubDownloader(self.token)
        self.org = self.cfg.HUB.ORG

    def _dev_model_path(self):
        pt = Paths()
        path = pt.models_development_path()
        if path is not None:
            path = os.path.join(path, self.model_id)
        if pt.exists(path):
            return path
        else:
            path = pt.ersilia_development_path()
            if path is not None:
                path = os.path.join(path, "test", "models", self.model_id)
            if pt.exists(path):
                return path
        return None

    @staticmethod
    def _copy_from_local(src, dst):
        shutil.copytree(src, dst)

    def _copy_from_github(self, dst):
        self.github_down.clone(org=self.org, repo=self.model_id, destination=dst)

    def get(self):
        """Copy model repository from local or download from GitHub"""
        folder = self._model_path(self.model_id)
        dev_model_path = self._dev_model_path()
        if dev_model_path is not None:
            self.logger.debug(
                "Copying from local {0} to {1}".format(dev_model_path, folder)
            )
            self._copy_from_local(dev_model_path, folder)
        else:
            self.logger.debug("Cloning from github to {0}".format(folder))
            self._copy_from_github(folder)


#  TODO: work outside GIT LFS
class ModelParametersGetter(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )

    @staticmethod
    def _requires_parameters(model_path):
        pf = PackFile(model_path)
        return pf.needs_model()

    def _get_destination(self):
        model_path = self._model_path(self.model_id)
        return os.path.join(model_path, MODEL_DIR)

    def get(self):
        """Create a ./model folder in the model repository"""
        model_path = self._model_path(self.model_id)
        folder = self._get_destination()
        if not os.path.exists(folder):
            os.mkdir(folder)
        if not self._requires_parameters(model_path):
            return None
        if not os.path.exists(folder):
            raise Exception


class ModelGetter(BaseAction):
    def __init__(self, model_id, repo_path, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.model_id = model_id
        self.repo_path = repo_path
        self.mrg = ModelRepositoryGetter(model_id=model_id, config_json=config_json)
        self.mpg = ModelParametersGetter(model_id=model_id, config_json=config_json)

    def _get_repository(self):
        self.mrg.get()

    def _get_model_parameters(self):
        self.mpg.get()

    def _copy_from_repo_path(self):
        dst = self._model_path(self.model_id)
        src = self.repo_path
        shutil.copytree(src, dst)

    def get(self):
        if self.repo_path is None:
            self._get_repository()
            self._get_model_parameters()
        else:
            self._copy_from_repo_path()
