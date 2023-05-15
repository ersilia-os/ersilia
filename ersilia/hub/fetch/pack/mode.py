import os
import json

from .... import ErsiliaBase
from ...bundle.repo import DockerfileFile
from ....utils.versioning import Versioner
from ....utils.system import SystemChecker
from ....setup.requirements.conda import CondaRequirement
from ....setup.requirements.docker import DockerRequirement
from ....default import MODEL_CONFIG_FILENAME

AVAILABLE_MODES = ["system", "venv", "conda", "docker"]


class PackModeDecision(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.versioner = Versioner(config_json=config_json)

    def _correct_protobuf(self, version, dockerfile, protobuf_version="3.19.5"):
        if version["version"] == "0.11.0":
            self.logger.debug(
                "Custom Ersilia BentoML is used, no need for modifying protobuf version"
            )
            return dockerfile
        if "0.11" in version["version"]:
            dockerfile.append_run_command(
                "pip install protobuf=={0}".format(protobuf_version)
            )
            self.logger.info(
                "Since BentoML is version 0.11, protobuf will been downgraded to {0}".format(
                    protobuf_version
                )
            )
        return dockerfile

    def decide_from_config_file_if_available(self):
        folder = self._model_path(self.model_id)
        if not os.path.exists(os.path.join(folder, MODEL_CONFIG_FILENAME)):
            return None
        with open(os.path.join(folder, MODEL_CONFIG_FILENAME), "r") as f:
            model_config = json.load(f)
        if "default_mode" in model_config:
            default_mode = model_config["default_mode"]
            if default_mode not in AVAILABLE_MODES:
                raise Exception(
                    "The model default_mode specified in the config.json file of the model repo is not correct. It should be one of {0}".format(
                        " ".join(AVAILABLE_MODES)
                    )
                )
            else:
                return default_mode
        return None

    def decide(self):
        sc = SystemChecker()
        if sc.is_github_action():
            self.logger.debug(
                "Code is being run inside a GitHub Actions workflow. Use conda as a by-default mode."
            )
            return "conda"
        mode = self.decide_from_config_file_if_available()
        if mode is not None:
            self.logger.debug("Mode is already specified in the model repository")
            self.logger.debug("Mode: {0}".format(mode))
            return mode
        folder = self._model_path(self.model_id)
        self.logger.debug(
            "Check if model can be run with vanilla (system) code (i.e. dockerfile has no installs)"
        )
        dockerfile = DockerfileFile(folder)
        self.logger.debug("Check bentoml and python version")
        version = dockerfile.get_bentoml_version()
        self.logger.info("BentoML version {0}".format(version))
        dockerfile = self._correct_protobuf(version, dockerfile)
        if not dockerfile.has_runs():
            same_python = version["python"] == self.versioner.python_version(
                py_format=True
            )
            same_bentoml = version["version"] == self.versioner.bentoml_version()
            if same_python and same_bentoml:
                self.logger.debug("Same python and same bentoml, run in system")
                self.logger.debug("Mode: system")
                return "system"
        self.logger.debug("Model needs some installs")
        cmds = dockerfile.get_install_commands()
        if cmds is None:
            self.logger.debug("No Dockerfile found...")
            raise Exception("No Dockerfile found!")
        self.logger.debug("Checking if only python/conda install will be sufficient")
        if cmds["exclusive_conda_and_pip"]:
            condareq = CondaRequirement()
            if not cmds["conda"] and not condareq.is_installed():
                self.logger.debug("Mode: venv")
                return "venv"
            else:
                self.logger.debug("Mode: conda")
                return "conda"
        else:
            self.logger.debug(
                "The python/conda installs may not be sufficient, trying docker"
            )
            self.logger.debug("Mode: docker")
            dockerreq = DockerRequirement()
            if dockerreq.is_inside_docker():
                return "conda"
            if dockerreq.is_installed():
                return "docker"
            else:
                return "conda"
