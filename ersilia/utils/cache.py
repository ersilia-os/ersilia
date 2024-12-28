import json
import threading

import redis.asyncio as redis

from ..default import REDIS_SERVER
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
            api_response = post_fn(missing_inputs)
            if not api_response:
                raise ValueError("API response is empty or invalid.")

            results_key = next(iter(api_response.keys()))
            metadata_key = "meta" if "meta" in api_response else None

            api_results = api_response[results_key]

            for input, api_result in zip(missing_inputs, api_results):
                key = f"{self.model_id}:{input['key']}"
                cached_value = {**input, **api_result}
                await self.client.setex(key, 3600, json.dumps(cached_value))  # 1h cache
                results.append(cached_value)

            metadata = api_response.get(metadata_key, {})
        else:
            metadata = {}

        return {"result": results, "meta": metadata}

    async def close(self):
        """
        Close the Redis client connection.
        """
        await self.client.aclose()
