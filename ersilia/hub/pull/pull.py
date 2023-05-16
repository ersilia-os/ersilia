from ... import ErsiliaBase
from ...utils.docker import SimpleDocker
from ...default import DOCKERHUB_ORG, DOCKERHUB_LATEST_TAG

import docker


class ModelPuller(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.simple_docker = SimpleDocker()
        self.model_id = model_id
        self.client = docker.from_env()

    def is_available_locally(self):
        return self.simple_docker.exists(
            DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG
        )

    def is_available(self):
        repository = DOCKERHUB_ORG + "/" + self.model_id
        tag = DOCKERHUB_LATEST_TAG
        try:
            self.client.images.get(f"{repository}:{tag}")
            return True
        except docker.errors.ImageNotFound:
            return False

    def pull(self):
        self.client.images.pull(
            "{0}/{1}".format(DOCKERHUB_ORG, self.model_id),
            tag=DOCKERHUB_LATEST_TAG,
            decode=True,
        )
