import os
import json
from re import I
from ... import ErsiliaBase
from ... import logger
from ...utils.paths import Paths
from ...utils.docker import SimpleDockerfileParser
from ...utils.conda import SimpleConda
from ...utils.system import SystemChecker
from ...default import (
    CONDA_ENV_YML_FILE,
    DOCKER_BENTO_PATH,
    DEFAULT_MODEL_ID,
    DOCKERFILE_FILE,
)

ROOT_CHECKFILE = "README.md"


class ReadmeFile(object):
    def __init__(self, path):
        self.path = os.path.abspath(path)

    def get_file(self):
        return os.path.join(self.path, "README.md")

    def check(self):
        return True


class ServiceFile(object):
    """checks the model service file

    Service File is needed to run a bentoml web app.
    Bentoml web app is named "service" by default.

    Attributes:
        path: directory path where model is stored. Uses os module.
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)

    def get_file(self):
        """gets the file from the specific directory"""
        return os.path.join(self.path, "src", "service.py")

    def _has_service_class(self):
        """checks if the service.py contains the Service Class."""
        search_string = "class Service("
        with open(self.get_file(), "r") as f:
            for l in f:
                if search_string == l[: len(search_string)]:
                    return True
        return False

    def rename_service(self):
        """renames the bentoml app from default "Service" to model_id."""
        ru = RepoUtils(self.path)
        model_id = ru.get_model_id()
        add_text = ru.rename_service(model_id)
        file_name = self.get_file()
        with open(file_name, "r") as f:
            text = f.read()
        if add_text in text:
            return
        text += os.linesep
        text += add_text
        with open(file_name, "w") as f:
            f.write(text)

    def add_info_api(self, information_file):
        """Adds and info api to the service"""
        file_name = self.get_file()
        with open(file_name, "r") as f:
            text = f.read()
        splitter_string = "Service.__name__"
        text = text.split(splitter_string)
        a = text[0]
        b = text[1]
        a += "    @api(input=JsonInput(), batch=True)\n"
        a += "    def info(self, input=None):\n"
        a += "        import json\n"
        a += "        data = json.load(open('{0}', 'r'))\n".format(information_file)
        a += "        return [data]\n\n"
        with open(file_name, "w") as f:
            s = a + splitter_string + b
            f.write(s)

    def check(self):
        """checks if the service.py contains the Service Class"""
        return self._has_service_class()


class PackFile(object):
    def __init__(self, path):
        self.path = os.path.abspath(path)

    def get_file(self):
        return os.path.join(self.path, "pack.py")

    def needs_model(self):
        # TODO: work on this function to account for more cases.
        file_name = self.get_file()
        line = None
        with open(file_name, "r") as f:
            for l in f:
                if ".pack(" in l:
                    line = l
        if line is None:
            return False
        if "None" in line:
            return False
        return True

    def check(self):
        return True


class DockerfileFile(object):
    """Checks the model Dockerfile.

    Dockerfile specifies the bentoml version,
    conda environment,
    python version.

    Attributes:
        path: directory path where model is stored. Uses os module.
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.parser = SimpleDockerfileParser(self.path)
        self.conda = SimpleConda()

    def get_file(self):
        """gets the file from the specific directory"""
        return os.path.join(self.path, DOCKERFILE_FILE)

    def get_bentoml_version(self):
        """Identifies the bentoml version required for the models.

        Uses the parser from dockerfile module to retrieve the baseimage

        Returns:
            A dictionary collecting bentoml version, slim (no conda env)
            and python version. For example:
            {"version":0.9.2, "slim" = True, "python":py37}
        """
        bimg = self.parser.baseimage
        bimg = bimg.split("/")
        if len(bimg) != 2:
            return None
        if bimg[0] != "bentoml":
            return None
        img = bimg[1]
        img = img.split(":")
        if len(img) != 2:
            return None
        if img[0] != "model-server":
            return None
        tag = img[1]
        if "-slim-" in tag:
            slim = True
        else:
            slim = False
        if slim:
            tag = tag.split("-")
            if len(tag) != 3:
                return None
        else:
            tag = tag.split("-")
            if len(tag) != 2:
                return None
        result = {"version": tag[0], "slim": slim, "python": tag[-1]}
        if result["python"] in [
            "py30",
            "py31",
            "py32",
            "py33",
            "py34",
            "py35",
            "py36",
            "py37",
            "py38",
            "py39",
        ]:
            if SystemChecker().is_arm64():
                logger.warning(
                    "This model is trying to install to a python version in ARM64 that is below 3.10. Changing to 3.10."
                )
                result["python"] = "py310"
        return result

    def has_runs(self):
        if self.parser.get_runs():
            return True
        else:
            return False

    def needs_conda(self):
        fn = self.get_file()
        if fn is None:
            return False
        with open(fn, "r") as f:
            for l in f:
                if "conda" in l:
                    return True
        return False

    def get_install_commands_from_dockerfile(self, fn):
        dp = SimpleDockerfileParser(fn)
        runs = dp.get_runs()
        return runs

    def get_install_commands(self):
        fn = self.get_file()
        if fn is None:
            return None
        needs_conda = self.needs_conda()
        runs = (
            self.conda.get_conda_and_pip_install_commands_from_dockerfile_if_exclusive(
                fn
            )
        )
        if runs:
            exclusive_conda_and_pip = True
        else:
            exclusive_conda_and_pip = False
            runs = self.get_install_commands_from_dockerfile(fn)
        result = {
            "conda": needs_conda,
            "commands": runs,
            "exclusive_conda_and_pip": exclusive_conda_and_pip,
        }
        return result

    def append_run_command(self, cmd):
        fn = self.get_file()
        if fn is None:
            return None
        R = []
        with open(fn, "r") as f:
            for l in f:
                R += [l]
        i_last = 0
        for i, r in enumerate(R):
            if r.startswith("RUN"):
                i_last = i
        j_add = 0
        for j, r in enumerate(R[(i_last + 1) :]):
            if r.startswith(" "):
                j_add = j + 1
            else:
                break
        i_last = i_last + j_add
        R = (
            R[: (i_last + 1)]
            + ["RUN {0}{1}".format(cmd, os.linesep)]
            + R[(i_last + 1) :]
        )
        with open(fn, "w") as f:
            for r in R:
                f.write(r)
        self.parser = SimpleDockerfileParser(self.path)

    def check(self):
        return True


class Integrity(object):
    def __init__(self, path):
        self.path = os.path.abspath(path)

    def _readme_file(self):
        return ReadmeFile(self.path).get_file()

    def _service_file(self):
        return ServiceFile(self.path).get_file()

    def _pack_file(self):
        return PackFile(self.path).get_file()

    def has_readme(self):
        if os.path.exists(self._readme_file()):
            return True
        else:
            return False

    def has_service(self):
        if os.path.exists(self._service_file()):
            return True
        else:
            return False

    def has_pack(self):
        if os.path.exists(self._pack_file()):
            return True
        else:
            return False


class RepoUtils(ErsiliaBase):
    def __init__(self, path, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.dockerhub_org = self.cfg.EXT.DOCKERHUB_ORG
        self.config_in_img = os.path.join(
            self.cfg.ENV.DOCKER.IMAGE_REPODIR, self.cfg.HUB.CONFIG_FILE
        )
        if os.path.isdir(path):
            self.path = os.path.normpath(os.path.abspath(path))
        else:
            self.path = os.path.normpath(os.path.dirname(os.path.abspath(path)))

    def _get_model_id_from_path(self):
        p = Paths()
        return p.model_id_from_path(self.path)

    def _get_model_id_from_config(self):
        if not os.path.exists(self.config_in_img):
            return None
        with open(self.config_in_img, "r") as f:
            cfg = json.load(f)
        return cfg["model_id"].replace("'", "").replace('"', "")

    def get_model_id(self):
        model_id = self._get_model_id_from_path()
        if model_id is None:
            model_id = self._get_model_id_from_config()
        if model_id is None:
            model_id = DEFAULT_MODEL_ID
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

    def _inside_docker(self):
        if os.path.exists(DOCKER_BENTO_PATH):
            return True
        else:
            return False

    def get_conda_env_yml_file(self):
        if self._inside_docker():
            return os.path.join(DOCKER_BENTO_PATH, CONDA_ENV_YML_FILE)
        else:
            root = self._root_path()
            if root is None:
                model_id = self.get_model_id()
                # try to find yml in bundles
                yml = os.path.join(
                    self._bundles_dir,
                    model_id,
                    self._get_latest_bundle_tag(model_id),
                    CONDA_ENV_YML_FILE,
                )
                if os.path.exists(yml):
                    return yml
                # try to find yml in bentoml
                yml = os.path.join(
                    self._bentoml_dir,
                    model_id,
                    self._get_latest_bentoml_tag(model_id),
                    CONDA_ENV_YML_FILE,
                )
                if os.path.exists(yml):
                    return yml
            else:
                return os.path.join(root, CONDA_ENV_YML_FILE)

    def get_docker_repo_image(self, model_id):
        return os.path.join(
            self.dockerhub_org, "{0}:{1}".format(model_id, self.cfg.ENV.DOCKER.REPO_TAG)
        )

    @staticmethod
    def rename_service(model_id):
        cmd = "Service.__name__ = '%s'\n" % model_id
        cmd += "%s = Service" % model_id
        return cmd
