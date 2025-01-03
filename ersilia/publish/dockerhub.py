from .. import ErsiliaBase
from ..db.environments.managers import DockerManager
from ..default import DOCKERHUB_LATEST_TAG, DOCKERHUB_ORG
from ..utils.terminal import run_command


class DockerHubUploader(ErsiliaBase):
    """
    Class for uploading Docker images to DockerHub.

    Parameters
    ----------
    model_id : str
        The ID of the model to be uploaded.
    config_json : str, optional
        Path to the configuration JSON file.

    Examples
    --------
    .. code-block:: python

        uploader = DockerHubUploader(
            model_id="model_id",
            config_json="path/to/config.json",
        )
        uploader.set_credentials(
            docker_user="username", docker_pwd="password"
        )
        uploader.upload()
    """

    def __init__(self, model_id: str, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.docker_user = None
        self.docker_pwd = None

    def set_credentials(self, docker_user: str, docker_pwd: str):
        """
        Set DockerHub credentials.

        Parameters
        ----------
        docker_user : str
            DockerHub username.
        docker_pwd : str
            DockerHub password.
        """
        self.docker_user = docker_user
        self.docker_pwd = docker_pwd

    def build_image(self):
        """
        Build the Docker image for the model.
        """
        dm = DockerManager()
        dm.build(self.model_id, self.docker_user, self.docker_pwd)

    def upload(self):
        """
        Upload the Docker image to DockerHub.
        """
        self.build_image()
        cmd = "docker login --password {0} --username {1}".format(
            self.docker_pwd, self.docker_user
        )
        run_command(cmd)
        cmd = "docker push {0}/{1}:{2}".format(
            DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG
        )
        run_command(cmd)
