import shutil
import subprocess
import time

import docker
import docker.errors
import psutil
import redis

from ..default import (
    DEFAULT_DOCKER_NETWORK_BRIDGE,
    DEFAULT_DOCKER_NETWORK_NAME,
    DEFAULT_REDIS_MEMORY_USAGE_FRACTION,
    REDIS_CONTAINER_NAME,
    REDIS_DATA_VOLUME,
    REDIS_HOST,
    REDIS_IMAGE,
    REDIS_PORT,
)
from ..utils.docker import set_docker_host
from ..utils.logging import logger


class SetupRedis:
    """
    A utility class to ensure a Redis server is running using Docker.

    This class uses the Docker SDK for Python to manage the Redis container.
    It checks for the availability of the Redis image, handles container lifecycle
    operations (restart, create) and ensures that a Docker volume is mounted for
    data persistence.

    Parameters
    ----------
    cache: bool
        Whether to use redis cache or not
    maxmemory: float
        Fraction of memory used by redis

    """

    def __init__(self, cache, maxmemory):
        self.cache = cache
        self.maxmemory = maxmemory
        try:
            if not self._is_amenable()[0]:
                return
            set_docker_host()
            self.client = docker.from_env()
            self.client.ping()
            logger.info("Docker is available and running.")
        except Exception:
            logger.error("Docker is not installed or running.")

        self.network = None
        self.volume = None

    def _get_max_memory_limit(self):
        total = psutil.virtual_memory().total

        if self.maxmemory:
            return int(total * self.maxmemory)
        return int(DEFAULT_REDIS_MEMORY_USAGE_FRACTION * total)

    def _configure_redis_memory_policy(self):
        client = redis.Redis(host=REDIS_HOST, port=REDIS_PORT)
        if client.ping():
            client.config_set("maxmemory", self._get_max_memory_limit())
            client.config_set("maxmemory-policy", "allkeys-lru")

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

        if not self._is_amenable()[0]:
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

        time.sleep(0.5)
        self._configure_redis_memory_policy()
        logger.info("Redis server setup completed.")

    def _is_amenable(self):
        if not self._check_docker_installed():
            logger.warning(
                "Docker is not installed in your system. Model result will not be cached!"
            )
            return (False, "Docker is not installed!")

        if not self._is_docker_active():
            logger.warning("Docker is not active. Model result will not be cached!")
            return (False, "Docker is not active!")

        if not self.cache:
            logger.warning("Caching is disabled. Model result will not be cached!")
            self._remove_container_if_exists()
            return (False, "Caching is disabled using flag!")

        if self.cache:
            logger.info("Caching is amenable")
            return (True, "Caching enabled!")
        logger.info("Caching is amenable")
        return (True, "Caching enabled!")

    def _is_docker_active(self):
        try:
            subprocess.run(
                ["docker", "info"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            return True
        except subprocess.CalledProcessError:
            return False
        except FileNotFoundError:
            return False

    def _check_docker_installed(self):
        if shutil.which("docker") is None:
            logger.warning(
                "Docker is not installed on this runner. Skipping Docker operations."
            )
            return False
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
            logger.info(f"Container {REDIS_CONTAINER_NAME} restarted successfully. ")
        except docker.errors.APIError:
            logger.error("Failed to restart container.")

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

    def remove_redis_image_if_exists(self):
        """
        Checks if the specified Docker image exists. If it does, remove it forcefully.
        """
        try:
            image = self.client.images.get(REDIS_IMAGE)
        except docker.errors.ImageNotFound:
            logger.warning(f"Image '{REDIS_IMAGE}' does not exist.")
            return
        except docker.errors.APIError as api_err:
            logger.warning(f"Error accessing Docker API: {api_err}")
            return
        except Exception as e:
            logger.warning(f"An unexpected error occurred: {e}")
            return

        try:
            self.client.images.remove(image=image.id, force=True)
            logger.debug(f"Image '{REDIS_IMAGE}' has been removed.")
        except docker.errors.APIError as api_err:
            logger.error(f"Failed to remove image '{REDIS_IMAGE}': {api_err}")
        except Exception as e:
            logger.error(f"An unexpected error occurred while removing the image: {e}")

    def _remove_container_if_exists(self):
        client = docker.from_env()
        try:
            container = client.containers.get(REDIS_CONTAINER_NAME)
        except docker.errors.NotFound:
            logger.error(f"Container '{REDIS_CONTAINER_NAME}' does not exist.")
            return
        except docker.errors.APIError as api_err:
            logger.error(f"Error accessing Docker API: {api_err}")
            return
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}")
            return

        try:
            container.remove(force=True)
            logger.error(f"Container '{REDIS_CONTAINER_NAME}' has been removed.")
        except docker.errors.APIError as api_err:
            logger.error(
                f"Failed to remove container '{REDIS_CONTAINER_NAME}': {api_err}"
            )
        except Exception as e:
            logger.error(
                f"An unexpected error occurred while removing the container: {e}"
            )

    def _start_new_container(self):
        self._create_docker_network()
        self._create_data_volume()

        try:
            container = self.client.containers.get(REDIS_CONTAINER_NAME)
            networks = container.attrs["NetworkSettings"]["Networks"].keys()
            if DEFAULT_DOCKER_NETWORK_NAME not in networks:
                logger.info(
                    f"Attaching existing Redis to network {DEFAULT_DOCKER_NETWORK_NAME}"
                )
                self.network.connect(container)
            else:
                logger.info(
                    "Redis container is already running on the correct network."
                )
            return container

        except docker.errors.NotFound:
            container = self.client.containers.run(
                image=REDIS_IMAGE,
                name=REDIS_CONTAINER_NAME,
                volumes={self.volume.name: {"bind": "/data", "mode": "rw"}},
                detach=True,
                ports={f"{REDIS_PORT}/tcp": REDIS_PORT},
                network=DEFAULT_DOCKER_NETWORK_NAME,
                auto_remove=True,
            )
            logger.info(
                f"Started new Redis container on {self.network.name}: {container.id}"
            )
            return container

    def _pull_and_start_container(self):
        try:
            logger.info(f"Pulling image {REDIS_IMAGE}...")
            self.client.images.pull(REDIS_IMAGE)
            self._start_new_container()
        except docker.errors.APIError:
            logger.error("Failed to pull image and start container.")
