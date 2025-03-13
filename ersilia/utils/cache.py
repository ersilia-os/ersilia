import json
import subprocess
import threading

import redis.asyncio as redis

from ..default import (
    DEFAULT_DOCKER_NETWORK_NAME,
    REDIS_CONTAINER_NAME,
    REDIS_IMAGE,
    REDIS_PORT,
    REDIS_SERVER,
)
from ..utils.logging import logger


class RedisClient:
    """
    Singleton class to manage Redis client with asynchronous connections using a connection pool.

    This class ensures that only one instance of the Redis client is created across the application.
    It supports appending a `model_id` as a prefix to Redis keys, enabling unique storage for
    computations associated with different models.

    Attributes
    ----------
    _instance : RedisClient
        The single instance of the RedisClient class.
    _lock : threading.Lock
        A lock to ensure thread-safe access to the singleton instance.
    pool : redis.ConnectionPool
        Redis connection pool that is used to manage Redis connections.
    client : redis.Redis
        Redis client object for interacting with Redis.
    model_id : str
        A string used as a prefix for all Redis keys to uniquely identify computations by model.
    """

    _instance = None
    _lock = threading.Lock()

    def __new__(cls):
        """
        Create a new instance of the class.

        Returns
        -------
        object
            The created instance.
        """
        if not cls._instance:
            with cls._lock:
                if not cls._instance:
                    cls._instance = super(RedisClient, cls).__new__(cls)
                    cls._instance.pool = redis.ConnectionPool.from_url(REDIS_SERVER)
                    cls._instance.client = redis.Redis(
                        connection_pool=cls._instance.pool
                    )
                    cls._instance.model_id = ""
        return cls._instance

    def set_model_id(self, model_id):
        """
        Set the model ID.

        Parameters
        ----------
        model_id : str
            The model ID to set.
        """
        self.model_id = model_id

    async def check_redis_server(self):
        """
        Check if the Redis server is running.

        Returns
        -------
        bool
            True if the Redis server is running, False otherwise.

        Raises
        ------
        ConnectionError
            If the Redis server is not reachable.
        """
        try:
            if await self.client.ping():
                logger.info("Redis server is running.")
                return True
        except Exception as e:
            raise ConnectionError("Redis server is not reachable.") from e

    async def get_conv_computed_result(self, input_data, post_fn, mapping, batch_size):
        """
        Get the computed result from Redis or compute it if not available.

        Parameters
        ----------
        input_data : list
            List of input data dictionaries.
        post_fn : function
            Function to call to compute the result if not available in Redis.

        Returns
        -------
        dict
            Dictionary containing the result and metadata.

        Raises
        ------
        ValueError
            If the model ID is not set or the API response is empty or invalid.
        """
        if not self.model_id:
            raise ValueError(
                "Model ID is not set. Use `set_model_id` to set a model ID before fetching results."
            )

        keys = [f"{self.model_id}:{input['key']}" for input in input_data]

        cached_results = await self.client.mget(keys)
        results = []
        missing_inputs = []

        for input, cached_result in zip(input_data, cached_results):
            if cached_result:
                logger.info(f"Cache hit for key: {input['key']}")
                cached_results_entry = json.loads(cached_result)
                logger.info(f"cached_results_entry: {cached_results_entry}")
                results.append(cached_results_entry)
            else:
                logger.info(f"Cache miss for key: {input['key']}")
                missing_inputs.append(input)

        if missing_inputs:
            logger.info(f"Missing inputs: {missing_inputs}")
            api_response = post_fn(batch_size, input_data, mapping)
            api_response = [
                {"value": item["output"]["value"]} for item in api_response.values()
            ]
            logger.info(f"API response: {api_response}")
            if not api_response:
                raise ValueError("API response is empty or invalid.")

            for input, api_result in zip(missing_inputs, api_response):
                key = f"{self.model_id}:{input['key']}"
                cached_value = {**input, **api_result}
                api_dict = {**api_result}
                await self.client.setex(
                    key, 604800, json.dumps(cached_value)
                )  # 1 week cache
                results.append(api_dict)

        results = {
            index: {
                "input": {
                    "key": item["key"],
                    "input": item["input"],
                    "text": item["text"],
                },
                "output": {"value": item["value"]},
            }
            for index, item in enumerate(results)
        }
        logger.info(f"Computed results: {results}")
        return results

    async def get_computed_result(self, input_data, post_fn):
        """
        Get the computed result from Redis or compute it if not available.

        Parameters
        ----------
        input_data : list
            List of input data dictionaries.
        post_fn : function
            Function to call to compute the result if not available in Redis.

        Returns
        -------
        dict
            Dictionary containing the result and metadata.

        Raises
        ------
        ValueError
            If the model ID is not set or the API response is empty or invalid.
        """
        if not self.model_id:
            raise ValueError(
                "Model ID is not set. Use `set_model_id` to set a model ID before fetching results."
            )

        keys = [f"{self.model_id}:{input['key']}" for input in input_data]

        cached_results = await self.client.mget(keys)
        results = []
        missing_inputs = []

        for input, cached_result in zip(input_data, cached_results):
            if cached_result:
                logger.info(f"Cache hit for key: {input['key']}")
                cached_results_entry = json.loads(cached_result)
                results.append(cached_results_entry)
            else:
                logger.info(f"Cache miss for key: {input['key']}")
                missing_inputs.append(input)

        if missing_inputs:
            logger.info(f"Missing inputs: {missing_inputs}")
            api_response = post_fn(missing_inputs)
            logger.info(f"API response: {api_response}")
            if not api_response:
                raise ValueError("API response is empty or invalid.")
            metadata_key = "meta" if "meta" in api_response else None

            for input, api_result in zip(missing_inputs, api_response):
                key = f"{self.model_id}:{input['key']}"
                cached_value = {**input, **api_result}
                logger.info(f"cached_value: {cached_value}")
                api_dict = {**api_result}
                await self.client.setex(
                    key, 604800, json.dumps(cached_value)
                )  # 1 week cache
                results.append(api_dict)

            metadata = api_response[0].get(metadata_key, {})
        else:
            metadata = {}
        results = {"result": results, "meta": metadata}
        logger.info(f"Computed results: {results}")
        return results

    async def close(self):
        """
        Close the Redis client connection.
        """
        await self.client.aclose()


class SetupRedis:
    """
    A utility class to ensure a Redis server is running using Docker.

    This class handles scenarios where the Redis server is not running
    by managing Docker containers for Redis. It checks for Docker installation,
    verifies the availability of Redis images, and manages container lifecycle
    operations such as restarting or creating a new container.

    Methods
    -------
    ensure_redis_running()
        Ensures that a Redis server is running by managing Docker containers.
    """

    def ensure_redis_running(self):
        """
        Ensure that a Redis server is running, starting a Docker container if necessary.

        This method handles the following scenarios:
        1. If Redis is running, do nothing.
        2. If the Redis image is available but the container is stopped, restart it.
        3. If no Redis container exists, create and start a new one.

        Raises
        ------
        RuntimeError
            If Docker is not installed or the container fails to start.
        """
        logger.info("Checking Redis server status...")

        try:
            self._check_docker_installed()

            if self._is_image_available():
                container_status = self._get_container_status()
                if container_status == "exited":
                    self._restart_container()
                elif not container_status:
                    self._start_new_container()
            else:
                self._pull_and_start_container()

            logger.info("Redis server setup completed.")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                "Failed to setup Redis container. Ensure Docker is installed and running."
            ) from e

    def _check_docker_installed(self):
        subprocess.run(["docker", "--version"], check=True, stdout=subprocess.DEVNULL)

    def _is_image_available(self):
        image_check = subprocess.run(
            ["docker", "images", "-q", REDIS_IMAGE],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )
        return bool(image_check.stdout.strip())

    def _get_container_status(self):
        container_check = subprocess.run(
            [
                "docker",
                "ps",
                "-a",
                "--filter",
                f"name={REDIS_CONTAINER_NAME}",
                "--format",
                "{{.State}}",
            ],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )
        return container_check.stdout.strip()

    def _restart_container(self):
        logger.info(f"Restarting stopped container: {REDIS_CONTAINER_NAME}")
        subprocess.run(["docker", "start", REDIS_CONTAINER_NAME], check=True)
        logger.info(f"Container {REDIS_CONTAINER_NAME} restarted successfully.")

    def _start_new_container(self):
        logger.info(f"Starting a new container: {REDIS_CONTAINER_NAME}")
        subprocess.run(
            [
                "docker",
                "run",
                "-d",
                "--name",
                REDIS_CONTAINER_NAME,
                "-p",
                f"{REDIS_PORT}:{REDIS_PORT}",
                "--network",
                DEFAULT_DOCKER_NETWORK_NAME,
                "--rm",
                REDIS_IMAGE,
            ],
            check=True,
        )
        logger.info(f"Container {REDIS_CONTAINER_NAME} started successfully.")

    def _pull_and_start_container(self):
        logger.info(f"Pulling image {REDIS_IMAGE}...")
        subprocess.run(["docker", "pull", REDIS_IMAGE], check=True)
        self._start_new_container()
