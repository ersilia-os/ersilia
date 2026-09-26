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

    This class provides a programmatic interface to fetch, serve, run and manage
    models from the Ersilia Model Hub. It mirrors the CLI commands of the same
    name and can be used as a context manager: entering the ``with`` block serves
    the model and leaving it closes the model.

    Parameters
    ----------
    model_id : str
        Identifier of the model, e.g. ``"eos4e40"``.
    verbose : bool, default=False
        Print logs to the terminal.

    Examples
    --------
    >>> from ersilia.api import Model
    >>> model = Model("eos4e40")
    >>> model.fetch()
    >>> with model:
    ...     df = model.run(["CCO", "c1ccccc1"])
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

        Downloads the model and its dependencies from the specified source. If no source
        is specified, defaults to DockerHub. If DockerHub is requested but Docker is not
        running, returns False instead of proceeding.

        Parameters
        ----------
        overwrite : bool, default=True
            Whether to overwrite existing local copy of the model.
        from_dir : str or None, default=None
            Local directory path to fetch model from. If specified, other sources are ignored.
        from_github : bool, default=False
            Fetch from GitHub.
        from_dockerhub : bool or None, default=None
            Fetch from DockerHub. Defaults to True when no source is specified.
        version : str, default="latest"
            Version/tag of the model to fetch.
        from_s3 : bool, default=False
            Fetch from Amazon S3.
        from_hosted : bool, default=False
            Fetch from a hosted URL.
        hosted_url : str or None, default=None
            URL for the hosted model service (used when from_hosted=True).
        verbose : bool or None, default=None
            Enable verbose logging.
        **kwargs
            Additional arguments passed through for forward compatibility.

        Returns
        -------
        bool or result
            False if DockerHub was requested but Docker is not running.
            Otherwise returns the result from the underlying fetch command.

        Examples
        --------
        From DockerHub (default), or from source on GitHub:

        >>> Model("eos4e40").fetch()
        >>> Model("eos4e40").fetch(from_github=True)
        """
        # infer default: if no source specified, use DockerHub
        if from_dockerhub is None:
            from_dockerhub = not any([from_dir, from_github, from_s3, from_hosted])

        if from_dockerhub and not self._is_docker_running():
            return False

        return fetch.fetch(
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

    def serve(
        self,
        port: int = None,
        track: bool = False,
        tracking_use_case: str = "local",
        enable_cache: bool = False,
        read_store: bool = False,
        write_store: bool = False,
        access: bool = None,
        nearest_neighbors: bool = False,
        max_cache_memory_frac: float = None,
        verbose_flag: bool = False,
    ):
        """
        Serve the model locally as an API, ready to receive predictions.

        Parameters
        ----------
        port : int, optional
            Port for the model server. If unspecified, a free port is chosen.
        track : bool, default=False
            Track runs (input/output stats, errors, timing) and send them to
            Ersilia's tracking bucket.
        tracking_use_case : str, default="local"
            Tracking use case when ``track`` is True. One of ``"local"``,
            ``"hosted"``, ``"self-service"`` or ``"test"``.
        enable_cache : bool, default=False
            Cache predictions in a local Redis container for 7 days.
        read_store : bool, default=False
            Read precalculated predictions from the Isaura store.
        write_store : bool, default=False
            Write predictions to the Isaura store.
        access : str, optional
            Visibility of predictions written to the Isaura store, ``"public"``
            or ``"private"``. Required with ``write_store``.
        nearest_neighbors : bool, default=False
            Use nearest-neighbour search when reading from the Isaura store.
        max_cache_memory_frac : float, optional
            Maximum fraction of system memory Redis may use. Recommended
            values are between 0.2 and 0.7.
        verbose_flag : bool, default=False
            Print logs to the terminal.

        Returns
        -------
        dict
            ``url`` (where the model is served), ``session`` and ``server``.

        Raises
        ------
        RuntimeError
            If the model is not found or ``max_cache_memory_frac`` is outside
            the recommended range.

        Examples
        --------
        >>> model = Model("eos4e40")
        >>> url = model.serve()["url"]
        >>> model.close()
        """
        self._url, self.session, self.SRV = serve.serve(
            self.model_id,
            port=port,
            track=track,
            tracking_use_case=tracking_use_case,
            enable_cache=enable_cache,
            read_store=read_store,
            write_store=write_store,
            access=access,
            nearest_neighbors=nearest_neighbors,
            max_cache_memory_frac=max_cache_memory_frac,
            verbose_flag=self.verbose_mode or verbose_flag,
        )
        return {
            "url": self._url,
            "session": self.session,
            "server": self.SRV,
        }

    def run(self, input_list, batch_size=1000):
        """
        Run the served model on a list of inputs.

        Parameters
        ----------
        input_list : list of str
            Inputs to the model, e.g. SMILES strings.
        batch_size : int, default=1000
            Number of inputs sent to the model server per batch.

        Returns
        -------
        pandas.DataFrame
            One row per input, with ``key`` and ``input`` columns followed by
            the model's output columns.

        Examples
        --------
        >>> with Model("eos4e40") as model:
        ...     df = model.run(["CCO", "c1ccccc1"])
        """
        return run.run(self.model_id, input_list, batch_size)

    def close(self):
        """
        Close the current session and clean up associated resources.

        Terminates the model server and removes the session file, freeing up system resources.

        Returns
        -------
        bool
            True if the session was successfully closed, False otherwise.
        """
        result = close.close(self.model_id)
        self._url = None
        return result

    def info(self):
        """
        Show information about the served model.

        Includes its identifiers, description, code and parameters links,
        Docker Hub image and supported architectures.

        Returns
        -------
        dict or None
            The model information.

        Raises
        ------
        RuntimeError
            If no model is served in the current session.
        """
        return info.info(self.model_id)

    def example(self, n_samples=5, mode="random"):
        """
        Generate example inputs for the served model.

        Parameters
        ----------
        n_samples : int, default=5
            Number of examples to generate. Ignored in ``"curated"`` mode.
        mode : str, default="random"
            ``"random"`` samples inputs at random, ``"deterministic"`` always
            returns the same inputs, and ``"curated"`` returns the model's own
            example file.

        Returns
        -------
        list or None
            The example inputs, or None if no model is served.

        Examples
        --------
        >>> with Model("eos4e40") as model:
        ...     inputs = model.example(n_samples=10)
        ...     df = model.run(inputs)
        """
        return example.example(n_samples, mode=mode)

    def delete(self):
        """
        Delete the model from local storage.

        Removes the model files and associated artifacts from the system.

        Returns
        -------
        bool
            True if the model was successfully deleted, False otherwise.

        Raises
        ------
        RuntimeError
            If the model cannot be deleted.
        """
        return delete.delete(self.model_id, verbose=self.verbose_mode)

    def is_fetched(self):
        """
        Check whether the model has been fetched on this machine.

        Returns
        -------
        bool
            True if the model is fetched, False otherwise.
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
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            echo(
                "❌ Docker is NOT running locally. Please start Docker to use Ersilia models.",
                fg="red",
            )
            return False

    def __enter__(self):
        self.serve()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class Catalog(object):
    """
    Browse the models available locally or in the Ersilia Model Hub.

    Unlike :class:`Model`, this class is not tied to a specific model. Its
    :meth:`catalog` method mirrors ``ersilia catalog`` but returns the result
    for programmatic use.

    Parameters
    ----------
    verbose : bool, default=False
        Print logs to the terminal.

    Examples
    --------
    >>> from ersilia.api import Catalog
    >>> df = Catalog().catalog(hub=True)
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
        List models, or show the model card of one model.

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
            The catalog, which is also printed to the terminal as a table (or
            as JSON when ``as_json`` is True).
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
