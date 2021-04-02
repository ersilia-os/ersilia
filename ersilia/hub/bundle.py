import os
import yaml
from ..core.base import ErsiliaBase
from ..default import CONDA_ENV_YML_FILE, DOCKERFILE_FILE
from .repo import DockerfileFile
from dockerfile_parse import DockerfileParser


class BundleEnvironmentFile(ErsiliaBase):

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.dir = os.path.abspath(self._get_bundle_location(model_id))
        self.path = os.path.join(self.dir, CONDA_ENV_YML_FILE)
        self.exists = os.path.exists(self.path)

    def get_file(self):
        return self.path

    def _is_not_pip(self, dep):
        if type(dep) is str:
            if dep != "pip":
                return True
        return False

    def needs_conda(self):
        if not self.exists:
            return False
        with open(self.path, "r") as f:
            data = yaml.safe_load(f)
            dependencies = data["dependencies"]
            for dep in dependencies:
                if self._is_not_pip(dep):
                    return True
                else:
                    return False
            return False

    def check(self):
        return True


class BundleDockerfileFile(ErsiliaBase):

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.dir = os.path.abspath(self._get_bundle_location(model_id))
        self.path = os.path.join(self.dir, DOCKERFILE_FILE)
        self.exists = os.path.exists(self.path)
        self.parser = DockerfileParser(path=self.path)

    def get_file(self):
        return self.path

    def get_bentoml_version(self):
        return DockerfileFile(path=self.path).get_bentoml_version()

    def set_to_slim(self):
        ver = self.get_bentoml_version()
        if not ver: return
        if ver["slim"]: return
        img = "bentoml/model-server:{0}-slim-{1}".format(ver["version"], ver["python"])
        self.parser.baseimage = img

    def check(self):
        return True
