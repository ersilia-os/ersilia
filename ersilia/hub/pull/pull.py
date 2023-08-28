import docker
import requests

from ... import ErsiliaBase
from ...utils.terminal import yes_no_input, run_command
from ... import throw_ersilia_exception
from ...utils.exceptions_utils.pull_exceptions import DockerImageNotAvailableError

from ...utils.docker import SimpleDocker
from ...default import DOCKERHUB_ORG, DOCKERHUB_LATEST_TAG


class ModelPuller(ErsiliaBase):
    def __init__(self, model_id, overwrite=None, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.simple_docker = SimpleDocker()
        self.model_id = model_id
        self.client = docker.from_env()
        self.image_name = "{0}/{1}:{2}".format(
            DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG
        )
        self.overwrite = overwrite

    def is_available_locally(self):
        is_available = self.simple_docker.exists(
            DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG
        )
        if is_available:
            self.logger.debug("Image {0} is available locally".format(self.image_name))
            return True
        else:
            self.logger.debug(
                "Image {0} is not available locally".format(self.image_name)
            )
            return False

    def is_available_in_dockerhub(self):
        url = "https://hub.docker.com/v2/repositories/{0}/{1}/tags/{2}".format(
            DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG
        )
        response = requests.get(url)
        if response.status_code == 200:
            self.logger.debug(
                "The docker image {0} exists in DockerHub".format(self.image_name)
            )
            return True
        else:
            self.logger.debug(
                "The docker image {0} does not exist in DockerHub".format(
                    self.image_name
                )
            )
            return False

    def _delete(self):
        self.logger.debug(
            "Deleting locally available image {0}".format(self.image_name)
        )
        self.simple_docker.delete(
            org=DOCKERHUB_ORG, img=self.model_id, tag=DOCKERHUB_LATEST_TAG
        )

    @throw_ersilia_exception
    def pull(self):
        if self.is_available_locally():
            if self.overwrite is None:
                do_pull = yes_no_input(
                    "Requested image {0} is available locally. Do you still want to fetch it? [Y/n]".format(
                        self.model_id
                    ),
                    default_answer="Y",
                )
            elif self.overwrite:
                do_pull = True
            else:
                do_pull = False
            if not do_pull:
                self.logger.info("Skipping pulling the image")
                return
            self._delete()
        else:
            self.logger.debug("Docker image of the model is not available locally")
        if self.is_available_in_dockerhub():
            self.logger.debug(
                "Pulling image {0} from DockerHub...".format(self.image_name)
            )
            try:
                img = self.client.images.get(
                    "{0}/{1}".format(DOCKERHUB_ORG, self.model_id)
                )
                self.logger.debug(
                    f"Size of image: {img.attrs['Size'] / (1024*1024):.2f} MB"
                )
            except:
                self.logger.warning("Could not obtain size of image")
            try:
                self.client.images.pull(
                    "{0}/{1}".format(DOCKERHUB_ORG, self.model_id),
                    tag=DOCKERHUB_LATEST_TAG,
                    decode=True,
                )
                self.logger.debug("Image pulled succesfully!")
            except:
                self.logger.warning(
                    "Conventional pull did not work, Ersilia is now forcing linux/amd64 architecture"
                )
                run_command(
                    "docker pull {0}/{1}:{2} --platform linux/amd64".format(
                        DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG
                    )
                )
            # except: #TODO add better error
            #    raise DockerImageArchitectureNotAvailableError(model=self.model_id)
        else:
            self.logger.info("Image {0} is not available".format(self.image_name))
            raise DockerImageNotAvailableError(model=self.model_id)
