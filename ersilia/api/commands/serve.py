from ... import ErsiliaModel, logger
from ...utils.session import register_model_session
from ..echo import echo


def serve(
    model: str,
    port: int = None,
    track: bool = False,
    tracking_use_case: str = "local",
    enable_cache: bool = True,
    read_store: bool = False,
    write_store: bool = False,
    access: bool = None,
    nearest_neighbors: bool = False,
    max_cache_memory_frac: float = None,
    verbose_flag: bool = False,
):
    """
    Serves a specified model as an API.

    Args
    -------
        model: The model ID to be served. Can either be the eos identifier or the slug identifier.
        port: The port to use when creating a model server. If unspecified, Ersilia looks for empty ports to use on the user's system.
        track: Whether the model's runs should be tracked to monitor model and system performance.
        tracking_use_case: If --track is true, this command allows specification of the tracking use case. Current options are: local, hosted, self-service and test.
        enable_cache: Toggle Redis-based local caching on or off. If enabled, the results from model APIs will be cached for 7 days.
        read_store: Specifies to read from isaura store
        write_store: Specifies to write from isaura store
        access: Specifies access level to write to isaura store
        nearest_neighbors: Specifies nearest neighbor search when reading from isaura store
        max_cache_memory_frac: Sets the maximum fraction of memory to use by Redis for caching. Recommended value 0.2-0.7.

    Returns
    -------
        Model ID, URL, SRV, Session, SRV, Session, Caching Mode Status, Local Cache Status, Tracking Status

    Raises
    -------
        RuntimeError: If the model/URL is not valid or not found,
        or if the maximum cache memory fraction is outside of the recommended range.

    """
    echo("Serving model. This process may take some time...", fg="blue")

    if verbose_flag:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)

    if max_cache_memory_frac is not None:
        if not (0.2 <= max_cache_memory_frac <= 0.7):
            raise RuntimeError(
                "Maximum fraction of memory to use by Redis for caching is outside of recommended range (0.2–0.7)."
            )

    mdl = ErsiliaModel(
        model,
        output_source=None,
        preferred_port=port,
        cache=enable_cache,
        maxmemory=max_cache_memory_frac,
        read_store=read_store,
        write_store=write_store,
        access=access,
        nearest_neighbors=nearest_neighbors,
    )

    if not mdl.is_valid():
        raise RuntimeError(f"Model {mdl.model_id} is not valid or not found.")

    track_runs = tracking_use_case if track else None

    mdl.serve(track_runs=track_runs)

    if mdl.url is None:
        raise RuntimeError("No URL found. Service unsuccessful.")

    register_model_session(mdl.model_id, mdl.session._session_dir)

    apis = mdl.get_apis()

    additional_apis = None
    if apis != ["run"]:
        additional_apis = []
        for api in apis:
            if api != "run":
                additional_apis.append(api)

    return mdl.url, mdl.session._session_dir, mdl.scl
