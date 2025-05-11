import json
import os

from ... import ErsiliaBase
from ...default import (
    CONDA_ENV_YML_FILE,
    DEFAULT_MODEL_ID,
    DOCKER_BENTO_PATH,
    DOCKERFILE_FILE,
)
from ...utils.conda import SimpleConda
from ...utils.docker import SimpleDockerfileParser
from ...utils.paths import Paths

ROOT_CHECKFILE = "README.md"


class ReadmeFile(object):
    """
    Class to handle README file operations.

    Parameters
    ----------
    path : str
        The directory path where the README file is located.
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)

    def get_file(self):
        """
        Get the path to the README file.

        Returns
        -------
        str
            The path to the README file.
        """
        return os.path.join(self.path, "README.md")

    def check(self):
        """
        Check if the README file exists.

        Returns
        -------
        bool
            True if the README file exists, False otherwise.
        """
        return True


class ServiceFile(object):
    """
    Class to handle service file operations for BentoML web apps.

    This class provides methods to get the path to the service file, rename the service,
    add an info API, and check if the service file contains the Service class.

    Parameters
    ----------
    path : str
        The directory path where the service file is located.
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)

    def get_file(self) -> str:
        """
        Get the path to the service file.

        Returns
        -------
        str
            The path to the service file.
        """
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
        """
        Rename the BentoML service from the default "Service" to the model ID.
        """
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

    def add_info_api(self, information_file: str):
        """
        Add an info API to the service.

        Parameters
        ----------
        information_file : str
            The path to the information file to be used by the info API.
        """
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

    def check(self) -> bool:
        """
        Check if the service file contains the Service class.

        Returns
        -------
        bool
            True if the service file contains the Service class, False otherwise.
        """
        return self._has_service_class()


class PackFile(object):
    """
    Class to handle pack file operations.

    This class provides methods to get the path to the pack file, check if the pack file needs a model,
    and check if the pack file exists.

    Parameters
    ----------
    path : str
        The directory path where the pack file is located.
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)

    def get_file(self) -> str:
        """
        Get the path to the pack file.

        Returns
        -------
        str
            The path to the pack file.
        """
        return os.path.join(self.path, "pack.py")

    def needs_model(self) -> bool:
        """
        Check if the pack file needs a model. Specifically this determines whether the "pack.py" file
        requires a model by checking if the file contains lines with the .pack() method and whether "None"
        is specified as an argument.

        Returns
        -------
        bool
        """
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
        """
        Check if the pack file exists.

        Returns
        -------
        bool
        """
        return True


class DockerfileFile(object):
    """
    Class to handle Dockerfile operations for models.

    This class provides methods to get the path to the Dockerfile, get the BentoML version required for the model,
    check if the Dockerfile contains RUN commands, check if the Dockerfile requires Conda, get installation commands,
    and more.

    Parameters
    ----------
    path : str
        The directory path where the Dockerfile is located.
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.parser = SimpleDockerfileParser(self.path)
        self.conda = SimpleConda()

    def get_file(self) -> str:
        """
        Get the path to the Dockerfile.

        Returns
        -------
        str
            The path to the Dockerfile.
        """
        return os.path.join(self.path, DOCKERFILE_FILE)

    def get_bentoml_version(self) -> dict:
        """
        Get the BentoML version required for the model.

        Returns
        -------
        dict
            A dictionary containing the BentoML version, slim flag, and Python version.
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
        # if result["python"] in [
        #     "py30",
        #     "py31",
        #     "py32",
        #     "py33",
        #     "py34",
        #     "py35",
        #     "py36",
        #     "py37",
        #     "py38",
        #     "py39",
        # ]:
        # if SystemChecker().is_arm64():
        #     logger.warning(
        #         "This model is trying to install to a python version in ARM64 that is below 3.10. Changing to 3.10."
        #     )
        #     result["python"] = "py310"
        return result

    def get_python_version(self) -> str:
        """
        Get the python version.

        Returns
        -------
        str
            Python version in the format py**
        """
        result = self.get_bentoml_version()
        return result["python"]

    def has_runs(self) -> bool:
        """
        Check if the Dockerfile contains RUN commands.

        Returns
        -------
        bool
            True if the Dockerfile contains RUN commands, False otherwise.
        """
        if self.parser.get_runs():
            return True
        else:
            return False

    def needs_conda(self) -> bool:
        """
        Check if the Dockerfile requires Conda.

        Returns
        -------
        bool
            True if the Dockerfile requires Conda, False otherwise.
        """
        fn = self.get_file()
        if fn is None:
            return False
        with open(fn, "r") as f:
            for l in f:
                if "conda" in l:
                    return True
        return False

    def get_install_commands_from_dockerfile(self, fn):
        """
        Get the install commands from the Dockerfile.

        Parameters
        ----------
        fn : str
            The path to the Dockerfile.

        Returns
        -------
        list
            The list of RUN commands from the Dockerfile.
        """
        dp = SimpleDockerfileParser(fn)
        runs = dp.get_runs()
        return runs

    def get_install_commands(self) -> dict:
        """
        Get the installation commands from the Dockerfile.

        Returns
        -------
        dict
            A dictionary containing Conda requirement, commands, and exclusive Conda and pip flag.
        """
        fn = self.get_file()
        if fn is None:
            return None
        needs_conda = self.needs_conda()
        runs = (
            self.conda.get_conda_and_pip_install_commands_from_dockerfile_if_exclusive(
                fn
            )
        )
        if runs is not None:
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

    def append_run_command(self, cmd: str):
        """
        Append a RUN command to the Dockerfile.

        Parameters
        ----------
        cmd : str
            The RUN command to append.
        """
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

    def check(self) -> bool:
        """
        Check if the Dockerfile exists.

        Returns
        -------
        bool
            True if the Dockerfile exists, False otherwise.
        """
        return True


class Integrity(object):
    """
    Class to check the integrity of the model files.

    It provides methods to check if the README file, service file, and pack file exist.

    Parameters
    ----------
    path : str
        The directory path where the model files are located.
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)

    def _readme_file(self):
        return ReadmeFile(self.path).get_file()

    def _service_file(self):
        return ServiceFile(self.path).get_file()

    def _pack_file(self):
        return PackFile(self.path).get_file()

    def has_readme(self) -> bool:
        """
        Check if the README file exists.

        Returns
        -------
        bool
            True if the README file exists, False otherwise.
        """
        if os.path.exists(self._readme_file()):
            return True
        else:
            return False

    def has_service(self) -> bool:
        """
        Check if the service file exists.

        Returns
        -------
        bool
            True if the service file exists, False otherwise.
        """
        if os.path.exists(self._service_file()):
            return True
        else:
            return False

    def has_pack(self) -> bool:
        """
        Check if the pack file exists.

        Returns
        -------
        bool
            True if the pack file exists, False otherwise.
        """
        if os.path.exists(self._pack_file()):
            return True
        else:
            return False


class RepoUtils(ErsiliaBase):
    """
    Utility class for handling repository operations.

    It provides methods to get the model ID, get the path to the Conda environment YAML file,
    get the Docker repository image, and rename the BentoML service.

    Parameters
    ----------
    path : str
        The directory path where the repository is located.
    config_json : dict, optional
        Configuration settings in JSON format.
    """

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

    def get_model_id(self) -> str:
        """
        Get the model ID from the repository.

        Returns
        -------
        str
            The model ID.
        """
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

    def get_conda_env_yml_file(self) -> str:
        """
        Get the path to the Conda environment YAML file.

        Returns
        -------
        str
            The path to the Conda environment YAML file.
        """
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

    def get_docker_repo_image(self, model_id: str) -> str:
        """
        Get the Docker repository image for the model.

        Parameters
        ----------
        model_id : str
            The ID of the model.

        Returns
        -------
        str
            The Docker repository image for the model.
        """
        return os.path.join(
            self.dockerhub_org, "{0}:{1}".format(model_id, self.cfg.ENV.DOCKER.REPO_TAG)
        )

    @staticmethod
    def rename_service(model_id: str) -> str:
        """
        Rename the BentoML service to the model ID.

        Parameters
        ----------
        model_id : str
            The ID of the model.

        Returns
        -------
        str
            The command to rename the service.
        """
        cmd = "Service.__name__ = '%s'\n" % model_id
        cmd += "%s = Service" % model_id
        return cmd
