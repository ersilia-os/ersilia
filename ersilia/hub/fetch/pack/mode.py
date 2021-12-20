from .... import ErsiliaBase
from ...bundle.repo import DockerfileFile
from ....utils.versioning import Versioner
from ....setup.requirements.conda import CondaRequirement
from ....setup.requirements.docker import DockerRequirement

AVAILABLE_MODES = ["system", "venv", "conda", "docker"]


class PackModeDecision(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.versioner = Versioner(config_json=config_json)

    def decide(self):
        folder = self._model_path(self.model_id)
        self.logger.debug(
            "Check if model can be run with vanilla (system) code (i.e. dockerfile has no installs)"
        )
        dockerfile = DockerfileFile(folder)
        self.logger.debug("Check bentoml and python version")
        version = dockerfile.get_bentoml_version()
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
        self.logger.debug("Checking if only python/conda install will be sufficient")
        if cmds is not None:
            condareq = CondaRequirement()
            if not cmds["conda"] and not condareq.is_installed():
                self.logger.debug("Mode: venv")
                return "venv"
            else:
                self.logger.debug("Mode: conda")
                return "conda"
        else:
            self.logger.debug(
                "The python/conda installs are not sufficient, use docker"
            )
            self.logger.debug("Mode: docker")
            dockerreq = DockerRequirement()
            assert dockerreq.is_installed()
            return "docker"
