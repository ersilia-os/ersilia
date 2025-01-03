import collections
import os

import yaml
from dockerfile_parse import DockerfileParser

from ...core.base import ErsiliaBase
from ...default import CONDA_ENV_YML_FILE, DOCKERFILE_FILE
from ...hub.fetch import MODEL_INSTALL_COMMANDS_FILE, REQUIREMENTS_TXT
from .repo import DockerfileFile


class BundleEnvironmentFile(ErsiliaBase):
    """
    Class to handle the environment file for a model bundle.

    Specifically provides methods to get the path, bundle Conda requirement,
    add model installation commands of environment file, and check if the environment file exists.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    config_json : dict, optional
        Configuration settings in JSON format.
    """

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.dir = os.path.abspath(self._get_bundle_location(model_id))
        self.path = os.path.join(self.dir, CONDA_ENV_YML_FILE)
        self.exists = os.path.exists(self.path)

    def get_file(self):
        """
        Get the path to the environment file.

        Returns
        -------
        str
            The path to the environment file.
        """
        return self.path

    def _is_not_pip(self, dep):
        if type(dep) is str:
            if dep != "pip":
                return True
        return False

    def needs_conda(self):
        """
        Check if the environment file requires Conda.

        Returns
        -------
        bool
            True if the environment file requires Conda, False otherwise.
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

    def add_model_install_commands(self):
        """
        Add model installation commands to the environment file.
        """
        f0 = self.path
        with open(f0, "r") as f:
            data = yaml.safe_load(f)
        f1 = os.path.join(
            self._get_bundle_location(self.model_id), MODEL_INSTALL_COMMANDS_FILE
        )
        with open(f1, "r") as f:
            for l in f:
                l = l.rstrip(os.linesep)
                if l[:13] == "conda install":
                    l = l[13:]
                    l = l.replace(" -y", "")
                    # find channel
                    if " -c " in l:
                        c = l.split(" -c ")[1].lstrip().split(" ")[0]
                    else:
                        c = None
                    # find dependency
                    if c is not None:
                        l = l.replace(" -c {0}".format(c), "")
                    d = l.strip()
                    if c not in data["channels"]:
                        data["channels"] += [c]
                    if d not in data["dependencies"]:
                        data["dependencies"] += [d]
        data_ = collections.OrderedDict()
        data_["name"] = data["name"]
        data_["channels"] = data["channels"]
        data_["dependencies"] = data["dependencies"]
        data = dict((k, v) for k, v in data_.items())
        with open(f0, "w") as f:
            yaml.safe_dump(data, f, sort_keys=False)

    def check(self):  # TODO: Removing this fucntion
        """
        Check if the environment file exists.

        Returns
        -------
        bool
        """
        return True


class BundleRequirementsFile(ErsiliaBase):
    """
    Class to handle the requirements file for a model bundle.

    Specifically provides methods to add model installation commands to the requirements file and check if the requirements file exists.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    config_json : dict, optional
        Configuration settings in JSON format.
    """

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.dir = os.path.abspath(self._get_bundle_location(model_id))
        self.path = os.path.join(self.dir, REQUIREMENTS_TXT)
        self.exists = os.path.exists(self.path)

    def add_model_install_commands(self):
        """
        Add model installation commands to the requirements file.
        """
        f0 = os.path.join(self._get_bundle_location(self.model_id), REQUIREMENTS_TXT)
        reqs = []
        with open(f0, "r") as f:
            for l in f:
                reqs += [l.strip(os.linesep)]
        f1 = os.path.join(
            self._get_bundle_location(self.model_id), MODEL_INSTALL_COMMANDS_FILE
        )
        with open(f1, "r") as f:
            for l in f:
                if "pip " in l:
                    r = l.rstrip(os.linesep).split(" ")[-1]
                    if r not in reqs:
                        reqs += [r]
        with open(f0, "w") as f:
            for l in reqs:
                f.write(l + os.linesep)

    def check(self):  # TODO: Removing this fucntion
        """
        Check if the requirements file exists.

        Returns
        -------
        bool
        """
        return True


class BundleDockerfileFile(ErsiliaBase):
    """
    Class to handle the Dockerfile for a model bundle.

    It specifically provides methods to get the path to the Dockerfile, get the BentoML version required for the model, and more.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    config_json : dict, optional
        Configuration settings in JSON format.
    """

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.dir = os.path.abspath(self._get_bundle_location(model_id))
        self.path = os.path.join(self.dir, DOCKERFILE_FILE)
        self.exists = os.path.exists(self.path)
        self.parser = DockerfileParser(path=self.path)

    def get_file(self) -> str:
        """
        Get the path to the Dockerfile.

        Returns
        -------
        str
            The path to the Dockerfile.
        """
        return self.path

    def get_bentoml_version(self) -> dict:
        """
        Get the BentoML version required for the model.

        Returns
        -------
        dict
            A dictionary containing the BentoML version, slim flag, and Python version.
        """
        return DockerfileFile(path=self.path).get_bentoml_version()

    def set_to_slim(self):
        """
        Set the Dockerfile to use the slim version of the BentoML image.
        """
        ver = self.get_bentoml_version()
        if not ver:
            return
        if ver["slim"]:
            return
        img = "bentoml/model-server:{0}-slim-{1}".format(ver["version"], ver["python"])
        self.parser.baseimage = img
        content = self.parser.content
        with open(self.path, "w") as f:
            f.write(content)

    def set_to_full(self):
        """
        Set the Dockerfile to use the full version of the BentoML image.
        """
        ver = self.get_bentoml_version()
        if not ver:
            return
        if ver["slim"]:
            img = "bentoml/model-server:{0}-{1}".format(ver["version"], ver["python"])
            self.parser.baseimage = img
        content = self.parser.content
        with open(self.path, "w") as f:
            f.write(content)

    def check(self):  # TODO: Removing this fucntion
        """
        Check if the Dockerfile exists.

        Returns
        -------
        bool
        """
        return True
