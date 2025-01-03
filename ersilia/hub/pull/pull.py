import asyncio
import os
import re
import subprocess

import requests

from ... import ErsiliaBase, throw_ersilia_exception
from ...default import DOCKERHUB_LATEST_TAG, DOCKERHUB_ORG
from ...utils.docker import SimpleDocker
from ...utils.exceptions_utils.pull_exceptions import (
    DockerConventionalPullError,
    DockerImageNotAvailableError,
)
from ...utils.logging import make_temp_dir
from ...utils.terminal import run_command, yes_no_input

PULL_IMAGE = os.environ.get("PULL_IMAGE", "Y")


class ModelPuller(ErsiliaBase):
    """
    ModelPuller is responsible for pulling models from DockerHub.

    Parameters
    ----------
    model_id : str
        The ID of the model to be pulled.
    overwrite : bool, optional
        Whether to overwrite existing files.
    config_json : dict, optional
        Configuration settings for the puller.

    Examples
    --------
    .. code-block:: python

        puller = ModelPuller(
            model_id="eosxxxx", config_json=config
        )
        await puller.async_pull()
    """

    def __init__(
        self,
        model_id: str,
        overwrite: bool = None,
        config_json: dict = None,
        docker_tag: str = None,
    ):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.simple_docker = SimpleDocker()
        self.model_id = model_id
        self.docker_tag = docker_tag or DOCKERHUB_LATEST_TAG
        self.image_name = "{0}/{1}:{2}".format(
            DOCKERHUB_ORG, self.model_id, self.docker_tag
        )
        self.overwrite = overwrite

    def is_available_locally(self) -> bool:
        """
        Check if the Docker image is available locally.

        Returns
        -------
        bool
            True if the image is available locally, False otherwise.
        """
        is_available = self.simple_docker.exists(
            DOCKERHUB_ORG, self.model_id, self.docker_tag
        )
        if is_available:
            self.logger.debug("Image {0} is available locally".format(self.image_name))
            return True
        else:
            self.logger.debug(
                "Image {0} is not available locally".format(self.image_name)
            )
            return False

    def is_available_in_dockerhub(self) -> bool:
        """
        Check if the Docker image is available in DockerHub.

        Returns
        -------
        bool
            True if the image is available in DockerHub, False otherwise.
        """
        url = "https://hub.docker.com/v2/repositories/{0}/{1}/tags/{2}".format(
            DOCKERHUB_ORG, self.model_id, self.docker_tag
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
            org=DOCKERHUB_ORG, img=self.model_id, tag=self.docker_tag
        )

    def _get_size_of_local_docker_image_in_mb(self) -> float:
        try:
            image_name = "{0}/{1}:{2}".format(
                DOCKERHUB_ORG, self.model_id, self.docker_tag
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

    @throw_ersilia_exception()
    async def async_pull(self):
        """
        Asynchronously pull the Docker image.
        """
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

                pull_command = (
                    f"docker pull {DOCKERHUB_ORG}/{self.model_id}:{self.docker_tag}"
                )

                process = await asyncio.create_subprocess_shell(
                    pull_command,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                )

                async def log_stream(stream, log_method):
                    async for line in stream:
                        log_method(line.decode().strip())

                await asyncio.gather(
                    log_stream(process.stdout, self.logger.info),
                    log_stream(process.stderr, self.logger.error),
                )

                await process.wait()

                if process.returncode != 0:
                    self.logger.error(
                        f"Pull command failed with return code {process.returncode}"
                    )
                    raise subprocess.CalledProcessError(
                        process.returncode, pull_command
                    )

                self.logger.debug("Image pulled successfully!")

            except DockerConventionalPullError:
                self.logger.warning(
                    "Conventional pull did not work, Ersilia is now forcing linux/amd64 architecture"
                )
                force_pull_command = f"docker pull {DOCKERHUB_ORG}/{self.model_id}:{self.docker_tag} --platform linux/amd64"

                process = await asyncio.create_subprocess_shell(
                    force_pull_command,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                )

                await asyncio.gather(
                    log_stream(process.stdout, self.logger.info),
                    log_stream(process.stderr, self.logger.error),
                )

                await process.wait()

                if process.returncode != 0:
                    self.logger.error(
                        f"Forced pull command failed with return code {process.returncode}"
                    )
                    raise subprocess.CalledProcessError(
                        process.returncode, force_pull_command
                    )

                self.logger.debug("Forced pull completed successfully!")

            size = self._get_size_of_local_docker_image_in_mb()
            if size:
                self.logger.debug("Size of image {0} MB".format(size))
            else:
                self.logger.warning("Could not obtain size of image")
            return size
        else:
            self.logger.info("Image {0} is not available".format(self.image_name))
            raise DockerImageNotAvailableError(model=self.model_id)

    @throw_ersilia_exception()
    def pull(self):
        """
        This method pulls the Docker image non-asynchronously.
        """
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
                        DOCKERHUB_ORG, self.model_id, self.docker_tag, tmp_file
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
                        DOCKERHUB_ORG, self.model_id, self.docker_tag
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
