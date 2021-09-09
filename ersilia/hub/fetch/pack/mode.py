from .... import ErsiliaBase
from ...bundle.repo import DockerfileFile
from ....utils.versioning import Versioner


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
        # model needs some installs
        cmds = dockerfile.get_install_commands()
        # only python/conda installs
        if cmds is not None:
            #  no conda is necessary
            if not cmds["conda"]:
                self.logger.debug("Mode: venv")
                return "venv"
            #  conda is necessary
            else:
                self.logger.debug("Mode: conda")
                return "conda"
        # python/conda installs are not sufficient, use dockerfile
        else:
            self.logger.debug("Mode: docker")
            return "docker"
