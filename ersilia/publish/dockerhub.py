from .. import ErsiliaBase
from ..db.environments.managers import DockerManager
from ..utils.terminal import run_command
from ..default import DOCKERHUB_ORG, DOCKERHUB_LATEST_TAG


class DockerHubUploader(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.docker_user = None
        self.docker_pwd = None

    def set_credentials(self, docker_user, docker_pwd):
        self.docker_user = docker_user
        self.docker_pwd = docker_pwd

    def build_image(self):
        dm = DockerManager()
        dm.build(self.model_id, self.docker_user, self.docker_pwd)

    def upload(self):
        self.build_image()
        cmd = "docker login --password {0} --username {1}".format(
            self.docker_pwd, self.docker_user
        )
        run_command(cmd)
        cmd = "docker push {0}/{1}:{2}".format(
            DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG
        )
        run_command(cmd)
