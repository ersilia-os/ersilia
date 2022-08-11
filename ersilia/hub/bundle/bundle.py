import os
import yaml
import collections
from ...core.base import ErsiliaBase
from ...default import CONDA_ENV_YML_FILE, DOCKERFILE_FILE
from ...hub.fetch import MODEL_INSTALL_COMMANDS_FILE, REQUIREMENTS_TXT
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

    def add_model_install_commands(self):
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

    def check(self):
        return True


class BundleRequirementsFile(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.dir = os.path.abspath(self._get_bundle_location(model_id))
        self.path = os.path.join(self.dir, REQUIREMENTS_TXT)
        self.exists = os.path.exists(self.path)

    def add_model_install_commands(self):
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
        ver = self.get_bentoml_version()
        if not ver:
            return
        if ver["slim"]:
            img = "bentoml/model-server:{0}-{1}".format(ver["version"], ver["python"])
            self.parser.baseimage = img
        content = self.parser.content
        with open(self.path, "w") as f:
            f.write(content)

    def check(self):
        return True
