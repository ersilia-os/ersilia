import os
import yaml
from ..core.base import ErsiliaBase
from ..default import CONDA_ENV_YML_FILE, DOCKERFILE_FILE
from .repo import DockerfileFile
from dockerfile_parse import DockerfileParser


class BundleEnvironmentFile(ErsiliaBase):
    """Analyses the environment.yml file created during model bundling

    environment.yml specifies the necessary dependencies; none,
    conda environment, pip dependencies

    Attributes:
        model_id: standard string id for Ersilia models, ie eos0xxx
        config_json:
        path: directory where environment.yml is stored
        exists: environment.yml has been created

    """
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.dir = os.path.abspath(self._get_bundle_location(model_id))
        self.path = os.path.join(self.dir, CONDA_ENV_YML_FILE)
        self.exists = os.path.exists(self.path)

    def get_file(self):
        """ Fetches environment.yml file"""
        return self.path

    def _is_not_pip(self, dep):
        """Checks that pip dependencies are not specified

        Args:
            dep:
        Returns:
            Boolean value:
                True if pip dependencies not required
                False if pip dependencies are required
        """
        if type(dep) is str:
            if dep != "pip":
                return True
        return False

    def needs_conda(self):
        """ Checks if conda environment is required

        Args:
            exists: environment.yml path

        Returns:
            Boolean value:
                True if conda env is required
                False if conda env is not required

        """
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
    """Analyses the Dockerfile created by bentoml during model fetching

    Attributes:
        model_id:standard string id for Ersilia models, ie eos0xxx
        config_json:
        path: directory where Dockerfile is stored
        exists: Dockerfile has been created
        parser: DockerfileParser function
    """
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.dir = os.path.abspath(self._get_bundle_location(model_id))
        self.path = os.path.join(self.dir, DOCKERFILE_FILE)
        self.exists = os.path.exists(self.path)
        self.parser = DockerfileParser(path=self.path)

    def get_file(self):
        """gets file from directory path"""
        return self.path

    def get_bentoml_version(self):
        """gets bentoml version necessary to run the model

        Uses the get_bentml_version function defined in the DockerfileFile class
        from repo.py

        Returns:
            bento_ml version as a string, for example:
            bentoml/model-server:0.9.2-py37
        """
        return DockerfileFile(path=self.path).get_bentoml_version()

    def set_to_slim(self):
        """modifies dockerfile to include the -slim- option in bentoml version

        -slim- option allows model deployment without conda dependencies

        Returns:
            bentoml version as a string, including the -slim-.
            For example: bentoml/model-server:0.9.2-slim-py37
        """
        ver = self.get_bentoml_version()
        if not ver: return
        if ver["slim"]: return
        img = "bentoml/model-server:{0}-slim-{1}".format(ver["version"], ver["python"])
        self.parser.baseimage = img

    def check(self):
        return True
