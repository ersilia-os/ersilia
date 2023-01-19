import sys
from . import BaseAction
from ...bundle.status import ModelStatus
from ....utils.terminal import run_command
from ....db.environments.localdb import EnvironmentDb
from ....setup.requirements.docker import DockerRequirement


class ModelToolizer(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.docker_org = self.cfg.EXT.DOCKERHUB_ORG
        self.model_status = ModelStatus()

    def pip_install(self, model_id):
        """Install the model and distribute as a python package"""
        self.logger.debug("Distributing as a python package with pip")
        bento = self._get_bundle_location(model_id)
        run_command([sys.executable, "-m", "pip", "install", bento])

    def dockerize(self, model_id):
        """Containerize model using bentoml with docker"""
        req = DockerRequirement()
        if not req.is_installed():
            self.logger.info("Cannot dockerize. Please make sure docker is installed.")
            return
        self.logger.debug("Dockerizing")
        tag = self.cfg.ENV.DOCKER.LATEST_TAG
        cmd = "bentoml containerize {1}:{2} -t {0}/{1}:{2}".format(
            self.docker_org, model_id, tag
        )
        run_command(cmd)
        # store docker in the local environment database
        db = EnvironmentDb(config_json=self.config_json)
        db.table = "docker"
        db.insert(
            model_id=model_id, env="{0}/{1}:{2}".format(self.docker_org, model_id, tag)
        )

    def toolize(self, do_pip, do_docker):
        self.logger.debug("Checking if model needs to be integrated to a tool")
        if do_pip:
            self.logger.debug("Integrating to pip")
            if not self.model_status.is_pip(self.model_id):
                self.pip_install(self.model_id)
        if do_docker:
            self.logger.debug("Integrating to docker")
            if not self.model_status.is_docker(self.model_id):
                self.dockerize(self.model_id)
