import shutil
import subprocess

import docker

from ..default import (
    DEFAULT_DOCKER_NETWORK_BRIDGE,
    DEFAULT_DOCKER_NETWORK_NAME,
    REDIS_CONTAINER_NAME,
    REDIS_DATA_VOLUME,
    REDIS_IMAGE,
    REDIS_PORT,
)
from ..utils.logging import logger


class SetupRedis:
    """
    A utility class to ensure a Redis server is running using Docker.

    This class uses the Docker SDK for Python to manage the Redis container.
    It checks for the availability of the Redis image, handles container lifecycle
    operations (restart, create) and ensures that a Docker volume is mounted for
    data persistence.
    """

    def __init__(self):
        try:
            if not self._check_docker_installed():
                return
            self.client = docker.from_env()
            self.client.ping()
            logger.info("Docker is available and running.")
        except Exception as e:
            logger.error("Docker is not installed or running.")
            raise RuntimeError("Docker is not installed or running.") from e

        self.network = None
        self.volume = None

    def ensure_redis_running(self):
        """
        Ensure that a Redis server is running by managing the Docker container.

        This method checks if the Redis image exists, whether a container with the
        specified name is available and running, and starts or restarts the container
        as needed.

        Raises:
            RuntimeError: If the container fails to start.
        """
        logger.info("Checking Redis server status...")
        if not self._check_docker_installed():
            return

        if self._is_image_available():
            container = self._get_container()
            if container:
                container.reload()
                if container.status == "exited" or container.status != "running":
                    self._restart_container(container)
                else:
                    logger.info("Redis container is already running.")
            else:
                self._start_new_container()
        else:
            self._pull_and_start_container()

        logger.info("Redis server setup completed.")

    def _check_docker_installed(self):
        if shutil.which("docker") is None:
            logger.warning(
                "Docker is not installed on this runner. Skipping Docker operations."
            )
            return False
        subprocess.run(["docker", "--version"], check=True, stdout=subprocess.DEVNULL)
        return True

    def _is_image_available(self):
        images = self.client.images.list(name=REDIS_IMAGE)
        available = len(images) > 0
        logger.info(f"Image {REDIS_IMAGE} available: {available}")
        return available

    def _get_container(self):
        containers = self.client.containers.list(
            all=True, filters={"name": REDIS_CONTAINER_NAME}
        )
        if containers:
            logger.info(f"Found existing container: {REDIS_CONTAINER_NAME}")
            return containers[0]
        return None

    def _restart_container(self, container):
        logger.info(f"Restarting container: {REDIS_CONTAINER_NAME}")
        try:
            container.start()
            container.reload()
            logger.info(
                f"Container {REDIS_CONTAINER_NAME} restarted successfully. Status: {container.status}"
            )
        except docker.errors.APIError as e:
            logger.error("Failed to restart container.")
            raise RuntimeError("Failed to restart Redis container.") from e

    def _create_docker_network(self):
        networks = self.client.networks.list(names=[DEFAULT_DOCKER_NETWORK_NAME])
        if networks:
            self.network = networks[0]
            logger.info(f"Docker network already exists: {self.network.name}")
        else:
            self.network = self.client.networks.create(
                name=DEFAULT_DOCKER_NETWORK_NAME,
                driver=DEFAULT_DOCKER_NETWORK_BRIDGE,
                attachable=True,
            )
            logger.info(f"Docker network created: {self.network.name}")

    def _create_data_volume(self):
        volumes = self.client.volumes.list(filters={"name": REDIS_DATA_VOLUME})
        if volumes:
            self.volume = volumes[0]
            logger.info(f"Volume already exists: {self.volume.name}")
        else:
            self.volume = self.client.volumes.create(name=REDIS_DATA_VOLUME)
            logger.info(f"Volume created: {self.volume.name}")

    def _start_new_container(self):
        self._create_docker_network()
        self._create_data_volume()

        logger.info(f"Starting a new container: {REDIS_CONTAINER_NAME}")
        try:
            container = self.client.containers.run(
                image=REDIS_IMAGE,
                name=REDIS_CONTAINER_NAME,
                detach=True,
                ports={f"{REDIS_PORT}/tcp": REDIS_PORT},
                network=DEFAULT_DOCKER_NETWORK_NAME,
                volumes={self.volume.name: {"bind": "/data", "mode": "rw"}},
                auto_remove=True,  # Container will be removed on stop, but volume persists.
            )
            logger.info(
                f"Container {REDIS_CONTAINER_NAME} started successfully. ID: {container.id}"
            )
        except docker.errors.APIError as e:
            logger.error("Failed to start new container.")
            raise RuntimeError("Failed to start new Redis container.") from e

    def _pull_and_start_container(self):
        try:
            logger.info(f"Pulling image {REDIS_IMAGE}...")
            self.client.images.pull(REDIS_IMAGE)
            self._start_new_container()
        except docker.errors.APIError as e:
            logger.error("Failed to pull image and start container.")
            raise RuntimeError("Failed to pull and start Redis container.") from e
