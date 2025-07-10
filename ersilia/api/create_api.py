from .commands import close, delete, example, fetch, info, run, serve


class ErsiliaAPI:
    """
    Python API wrapper for interacting with Ersilia Model Hub models.

    This class provides a programmatic interface to run, serve, fetch, and manage
    machine learning models from the Ersilia Model Hub. It wraps the existing CLI-based
    functionality in a clean Pythonic interface and supports use in context managers.

    Parameters
    ----------
    model_id : str
        The unique identifier of the model to be managed and executed.

    Methods
    -------
    fetch():
        Downloads the specified model and its dependencies from DockerHub.

    serve():
        Serves the specified model locally to prepare for inference.

    run(input, output, batch_size):
        Runs the model on the given input and writes predictions to the output file.

    close():
        Terminates the model server and cleans up associated resources.

    info():
        Prints metadata and technical details about the specified model.

    example(file_name, simple, random, n_samples, deterministic):
        Generates example input files for the model.

    delete():
        Deletes the specified model and its local artifacts.

    __enter__():
        Context manager entry — automatically serves the model.

    __exit__(exc_type, exc_value, traceback):
        Context manager exit — automatically closes the served model.

    Usage Example
    -------------
    >>> from ersilia.api.create_api import ErsiliaAPI
    >>> molecular_weight = ErsiliaAPI("eos3b5e")
    >>> molecular_weight.fetch()
    >>> molecular_weight.serve()
    >>> with molecular_weight as model:
    >>>     model.info()
    """
    def __init__(self, model_id):
        self.model_id = model_id

    def fetch(self):
        """
        Fetches an Ersilia model to run it locally.

        This command allows users to fetch a specified model from the model hub (dockerhub, repo, s3 etc...).

        Returns
        -------
        Function: The fetch command function to be used by the API.
        Str: Confirmation message on success or warning message on failure.

        Raises
        -------
        RuntimeError: If both BentoML and FastAPI are used together.

        """
        fetch.fetch(
            model=self.model_id,
            overwrite=True,
            from_dir=None,
            from_github=False,
            from_dockerhub=True,
            version="latest",
            from_s3=False,
            from_hosted=False,
            with_fastapi=True,
            with_bentoml=None,
            hosted_url=None,
            verbose=False,
        )

    def serve(self):
        """
        Serves a specified model as an API.

        Args
        -------
            model: The model ID to be served. Can either be the eos identifier or the slug identifier.
            port: The port to use when creating a model server. If unspecified, Ersilia looks for empty ports to use on the user's system.
            track: Whether the model's runs should be tracked to monitor model and system performance.
            tracking_use_case: If --track is true, this command allows specification of the tracking use case. Current options are: local, hosted, self-service and test.
            enable_local_cache: Toggle Redis-based local caching on or off. If enabled, the results from model APIs will be cached for 7 days.
            local_cache_only: Specifies to fetch stored model results from local cache. The local caching system is powered by Redis.
            cloud_cache_only: Specifies to fetch stored model results from cloud cache. This allows to fetch model precalculated results in csv file in Ersilia model output format.
            cache_only: Specifies to fetch stored model results from both local and cloud cache. More details are given in a dump CLI.
            max_cache_memory_frac: Sets the maximum fraction of memory to use by Redis for caching. Recommended value 0.2-0.7.

        Returns
        -------
            Model ID, URL, SRV, Session, SRV, Session, Caching Mode Status, Local Cache Status, Tracking Status

        Raises
        -------
            RuntimeError: If the model/URL is not valid or not found,
            or if the maximum cache memory fraction is outside of the recommended range.

        """
        serve.serve(
            self.model_id,
            port=None,
            track=False,
            tracking_use_case="local",
            enable_local_cache=True,
            local_cache_only=False,
            cloud_cache_only=False,
            cache_only=False,
            max_cache_memory_frac=None,
            verbose=False,
        )

    def run(self, input, output, batch_size):
        """
        Runs the current model on a list of SMILES strings and
        returns the prediction as a pandas data frame.

        Args
        ----
        input: a list or a path to a CSV file containing SMILES strings.
        output: path to the output file where predictions will be saved.
        batch_size: number of SMILES to process per batch
        
        Returns
        -------
        function
            The run command function to be used by the API.
            A pandas df with the predictions.

        """
        print(run.run(self.model_id, input, output, batch_size))

    def close(self):
        close.close(self.model_id)

    def info(self):
        info.info(self.model_id)

    def example(self, file_name, simple, random, n_samples, deterministic):
        print(
            example.example(
                self.model_id, file_name, simple, random, n_samples, deterministic
            )
        )

    def delete(self):
        delete.delete(self.model_id, verbose=False)

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
