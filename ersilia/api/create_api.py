import subprocess

from .commands import (
    catalog,
    close,
    delete,
    example,
    fetch,
    info,
    is_fetched,
    run,
    serve,
)
from .echo import echo


class Model(object):
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
    >>> from ersilia.api import Model
    >>> molecular_weight = Model("eos3b5e")
    >>> molecular_weight.fetch()
    >>> molecular_weight.serve()
    >>> with molecular_weight as model:
    >>>     model.info()
    """

    def __init__(self, model_id, verbose=False):
        self.model_id = model_id
        self.verbose_mode = verbose
        self._url = None
        self.session = None
        self.SRV = None

    def fetch(
        self,
        *,
        overwrite: bool = True,
        from_dir: str | None = None,
        from_github: bool = False,
        from_dockerhub: bool | None = None,
        version: str = "latest",
        from_s3: bool = False,
        from_hosted: bool = False,
        hosted_url: str | None = None,
        verbose: bool | None = None,
        **kwargs,
    ):
        """
        Fetch an Ersilia model to run it locally.

        Parameters mirror the CLI:
        - overwrite: overwrite existing local copy (default: True)
        - from_dir: fetch from a local directory path
        - from_github: fetch from GitHub
        - from_dockerhub: fetch from DockerHub (if None, inferred)
        - version: version tag (default: "latest")
        - from_s3: fetch from S3
        - from_hosted: fetch from a hosted URL
        - hosted_url: URL to hosted assets (used when from_hosted=True)
        - verbose: override for verbose mode
        - **kwargs: passthrough for forward compatibility
        """
        # infer default: if no source specified, use DockerHub
        if from_dockerhub is None:
            from_dockerhub = not any([from_dir, from_github, from_s3, from_hosted])

        fetch.fetch(
            model=self.model_id,
            overwrite=overwrite,
            from_dir=from_dir,
            from_github=from_github,
            from_dockerhub=from_dockerhub,
            version=version,
            from_s3=from_s3,
            from_hosted=from_hosted,
            hosted_url=hosted_url,
            verbose_flag=(self.verbose_mode if verbose is None else verbose),
            **kwargs,
        )

    def serve(self, verbose=None):
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
        self._url, self.session, self.SRV = serve.serve(
            self.model_id,
            port=None,
            track=False,
            tracking_use_case="local",
            enable_local_cache=True,
            local_cache_only=False,
            cloud_cache_only=False,
            cache_only=False,
            max_cache_memory_frac=None,
            verbose_flag=self.verbose_mode or verbose,
        )

    def run(self, input_list, batch_size=1000):
        """
        Runs the current model on a list of input strings and
        returns the prediction as a pandas dataframe.

        Args
        ----
        input_list: a list containing input strings.
        batch_size: number of input strings to process per batch

        Returns
        -------
        function
            The run command function to be used by the API.
            A pandas df with the predictions.

        """
        return run.run(self.model_id, input_list, batch_size)

    def close(self):
        """
        This command closes the current session of the served model and cleans up any associated resources.

        Args
        -------
            model_id (str): ID of the model to delete.

        Returns
        -------
            str: Confirmation message on success or warning message on failure.

        Raises
        -------
            RuntimeError: If no model was served in the current session.
        """
        close.close(self.model_id)
        self._url = None

    def info(self):
        """
        Provides information about a specified model.

        This command allows users to get detailed information about a current active session,
        including information about Model Identifiers, Code and Parameters, Docker Hub link and Architectures.

        Args
        -------
        model_id (str): ID of the model to delete.

        Returns
        -------
        function: The info command function to be used by the API.
        str: Confirmation message on success or warning message on failure.

        Raises
        -------
        RuntimeError: If no model was served in the current session.
        """
        return info.info(self.model_id)

    def example(self, n_samples=5, mode="random"):
        """
        This command can sample inputs for a given model.

        Args
        -------
        model: The model ID to be served. Can either be the eos identifier or the slug identifier.
        simple: Simple inputs only contain the input column, while complete inputs also include key and the input.
        random: If the model source contains an example input file, when the predefined flag is set, then inputs are sampled from that file. Only the number of samples present in the file are returned, especially if --n_samples is greater than that number. By default, Ersilia samples inputs randomly.
        n_samples: Specify the number of example inputs to generate for the given model.
        deterministic: Used to generate examples data deterministically instead of random sampling. This allows when every time you run with example command with this flag you get the same types of examples.

        Returns
        -------
        Function: The exmaple command function to be used by the API.
        Str: Error message if no model was served in the current session.

        """
        return example.example(n_samples, mode=mode)

    def delete(self):
        """
        Deletes a specified model from local storage.

        Args:
            model_id (str): ID of the model to delete.

        Returns:
            str: Confirmation message on success or warning message on failure.

        Raises:
            RuntimeError: If the model cannot be deleted.
        """
        delete.delete(self.model_id, verbose=self.verbose_mode)

    def is_fetched(self):
        """
        Checks whether the model has been successfully fetched.

        Returns
        -------
        Echo Message indicating fetch status.
        """
        return is_fetched.is_fetched(self.model_id)

    def _is_docker_running(self):
        """
        Checks if Docker is running locally by calling `docker info`.

        Returns
        -------
        bool
            True if Docker is running, False otherwise.
        """
        try:
            subprocess.run(
                ["docker", "info"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )
            echo("✅ Docker is running locally.", fg="green")
        except (subprocess.CalledProcessError, FileNotFoundError):
            echo(
                "❌ Docker is NOT running locally. Please start Docker to use Ersilia models.",
                fg="red",
            )

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class Catalog(object):
    """
    This class enables users to browse and retrieve information about all models
    available in the Ersilia Hub. It is designed for general catalog-level operations
    and is not tied to any specific model instance. Use this class when you want to
    explore the hub, list available models,

    Typical usage includes listing the model catalog via the `catalog()` method,
    which mirrors the CLI behavior but returns a structured DataFrame for programmatic use.

    """

    def __init__(self, verbose=False):
        self.verbose_mode = verbose

    def catalog(
        self,
        hub=False,
        file_name=None,
        more=False,
        card=False,
        model=None,
        as_json=False,
        verbose=False,
    ):
        """
        API-compatible version of the catalog command with echo-based output.
        Parameters
        ----------
        hub : bool, default=False
            If True, fetch the catalog from the hub.
            If False, fetch the catalog from the local directory.
        file_name : str or None, default=None
            If specified, write the catalog to this file.
        more : bool, default=False
            If True, show more detail in catalog.
        card : bool, default=False
            If True, display the model card for a given model.
        model : str or None, default=None
            The model ID for which to display metadata.
        as_json : bool, default=False
            If True, return JSON output instead of a formatted table.
        verbose : bool, default=False
            If True, enable verbose logging.
        Returns
        -------
        pandas.DataFrame or dict or None
            A DataFrame containing the last two columns of the model catalog.
            Also prints the full catalog as a table or JSON to the terminal, depending on `as_json`.
        """
        return catalog.catalog(
            hub=hub,
            file_name=file_name,
            more=more,
            card=card,
            model=model,
            as_json=as_json,
            verbose=verbose,
        )
