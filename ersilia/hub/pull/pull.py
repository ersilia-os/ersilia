import requests
import subprocess
import tempfile
import json
import os
import re

from ... import ErsiliaBase
from ...utils.terminal import yes_no_input, run_command
from ... import throw_ersilia_exception
from ...utils.exceptions_utils.pull_exceptions import (
    DockerImageNotAvailableError,
    DockerConventionalPullError,
)

from ...utils.docker import SimpleDocker
from ...default import DOCKERHUB_ORG, DOCKERHUB_LATEST_TAG, EOS, MODEL_SIZE_FILE
from ...utils.logging import make_temp_dir

PULL_IMAGE = os.environ.get("PULL_IMAGE", "Y")


class ModelPuller(ErsiliaBase):
    def __init__(self, model_id, overwrite=None, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.simple_docker = SimpleDocker()
        self.model_id = model_id
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

    def _get_size_of_local_docker_image_in_mb(self):
        try:
            image_name = "{0}/{1}:{2}".format(
                DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG
            )
            result = subprocess.check_output(
                ["docker", "image", "inspect", image_name, "--format", "{{.Size}}"]
            )
            # Convert bytes to MB
            size_in_mb = int(result.strip()) / (1024 * 1024)
            return size_in_mb
        except subprocess.CalledProcessError:
            self.logger.warning("Image not found locally")
            return None

    @throw_ersilia_exception
    def pull(self):
        if self.is_available_locally():
            if self.overwrite is None:
                do_pull = yes_no_input(
                    "Requested image {0} is available locally. Do you still want to fetch it? [Y/n]".format(
                        self.model_id
                    ),
                    default_answer=PULL_IMAGE,
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
                self.logger.debug(
                    "Trying to pull image {0}/{1}".format(DOCKERHUB_ORG, self.model_id)
                )
                tmp_file = os.path.join(
                    make_temp_dir(prefix="ersilia-"), "docker_pull.log"
                )
                self.logger.debug("Keeping logs of pull in {0}".format(tmp_file))
                run_command(
                    "docker pull {0}/{1}:{2} > {3} 2>&1".format(
                        DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG, tmp_file
                    )
                )
                with open(tmp_file, "r") as f:
                    pull_log = f.read()
                    self.logger.debug(pull_log)
                if re.search(r"no match.*manifest", pull_log):
                    self.logger.warning(
                        "No matching manifest for image {0}".format(self.model_id)
                    )
                    raise DockerConventionalPullError(model=self.model_id)
                self.logger.debug("Image pulled succesfully!")
            except DockerConventionalPullError:
                self.logger.warning(
                    "Conventional pull did not work, Ersilia is now forcing linux/amd64 architecture"
                )
                run_command(
                    "docker pull {0}/{1}:{2} --platform linux/amd64".format(
                        DOCKERHUB_ORG, self.model_id, DOCKERHUB_LATEST_TAG
                    )
                )
            size = self._get_size_of_local_docker_image_in_mb()
            if size:
                self.logger.debug("Size of image {0} MB".format(size))
                # path = os.path.join(self._model_path(self.model_id), MODEL_SIZE_FILE)
                # with open(path, "w") as f:
                #     json.dump({"size": size, "units": "MB"}, f, indent=4)
                # self.logger.debug("Size written to {}".format(path))
            else:
                self.logger.warning("Could not obtain size of image")
            return size
        else:
            self.logger.info("Image {0} is not available".format(self.image_name))
            raise DockerImageNotAvailableError(model=self.model_id)
