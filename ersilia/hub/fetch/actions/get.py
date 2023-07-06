import os
import shutil
import tempfile
import zipfile
from . import BaseAction
from ....utils.download import GitHubDownloader, S3Downloader
from ....utils.paths import Paths
from ...bundle.repo import PackFile, DockerfileFile
from ....utils.exceptions_utils.throw_ersilia_exception import throw_ersilia_exception
from ....utils.exceptions_utils.fetch_exceptions import (
    FolderNotFoundError,
    S3DownloaderError,
)

from ....default import S3_BUCKET_URL_ZIP

MODEL_DIR = "model"


class ModelRepositoryGetter(BaseAction):
    def __init__(self, model_id, config_json, force_from_github, force_from_s3):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.token = self.cfg.HUB.TOKEN
        self.github_down = GitHubDownloader(self.token)
        self.s3_down = S3Downloader()
        self.org = self.cfg.HUB.ORG
        self.force_from_github = force_from_github
        self.force_from_s3 = force_from_s3

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

    def _copy_zip_from_s3(self, dst):
        self.logger.debug("Downloading model from S3 in zipped format")
        tmp_file = os.path.join(tempfile.mkdtemp("ersilia-"), "model.zip")
        self.s3_down.download_from_s3(
            bucket_url=S3_BUCKET_URL_ZIP,
            file_name=self.model_id + ".zip",
            destination=tmp_file,
        )
        self.logger.debug("Extracting model from {0}".format(tmp_file))
        dst = "/".join(dst.split("/")[:-1])
        self.logger.debug("...to {0}".format(dst))
        with zipfile.ZipFile(tmp_file, "r") as zip_ref:
            zip_ref.extractall(dst)

    def _change_py_version_in_dockerfile_if_necessary(self):
        path = self._model_path(model_id=self.model_id)
        df = DockerfileFile(path=path)
        version = df.get_bentoml_version()
        self.logger.debug(version)
        dockerfile_path = os.path.join(path, "Dockerfile")
        with open(dockerfile_path, "r") as f:
            R = f.readlines()
        S = []
        for r in R:
            if r.startswith("FROM "):
                r = r.split("-")
                if r[-1].startswith("py"):
                    p = version["python"]
                    r = "-".join(r[:-1] + [p])
                else:
                    r = "-".join(r)
            S += [r]
        with open(dockerfile_path, "w") as f:
            for s in S:
                f.write(s + os.linesep)

    @throw_ersilia_exception
    def get(self):
        """Copy model repository from local or download from S3 or GitHub"""
        folder = self._model_path(self.model_id)
        dev_model_path = self._dev_model_path()
        if dev_model_path is not None:
            self.logger.debug(
                "Copying from local {0} to {1}".format(dev_model_path, folder)
            )
            self._copy_from_local(dev_model_path, folder)
        else:
            if self.force_from_github:
                self._copy_from_github(folder)
            else:
                try:
                    self.logger.debug("Trying to download from S3")
                    self._copy_zip_from_s3(folder)
                except:
                    self.logger.debug(
                        "Could not download in zip format in S3. Downloading from GitHub repository."
                    )
                    if self.force_from_s3:
                        raise S3DownloaderError(model_id=self.model_id)
                    else:
                        self._copy_from_github(folder)
        self._change_py_version_in_dockerfile_if_necessary()


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

    @throw_ersilia_exception
    def get(self):
        """Create a ./model folder in the model repository"""
        model_path = self._model_path(self.model_id)
        folder = self._get_destination()
        if not os.path.exists(folder):
            os.mkdir(folder)
        if not self._requires_parameters(model_path):
            return None
        if not os.path.exists(folder):
            raise FolderNotFoundError(folder)


class ModelGetter(BaseAction):
    def __init__(
        self, model_id, repo_path, config_json, force_from_gihtub, force_from_s3
    ):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.model_id = model_id
        self.repo_path = repo_path
        self.mrg = ModelRepositoryGetter(
            model_id=model_id,
            config_json=config_json,
            force_from_github=force_from_gihtub,
            force_from_s3=force_from_s3,
        )
        self.mpg = ModelParametersGetter(model_id=model_id, config_json=config_json)

    def _get_repository(self):
        self.mrg.get()

    def _get_model_parameters(self):
        self.mpg.get()

    def _copy_from_repo_path(self):
        dst = self._model_path(self.model_id)
        src = self.repo_path
        shutil.copytree(src, dst)

    @throw_ersilia_exception
    def get(self):
        if self.repo_path is None:
            self._get_repository()
            self._get_model_parameters()
        else:
            self._copy_from_repo_path()
        if not os.path.exists(self._model_path(self.model_id)):
            raise FolderNotFoundError(os.path.exists(self._model_path(self.model_id)))
