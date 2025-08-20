from ... import ErsiliaModel, logger
from ...store.utils import OutputSource
from ...utils.cache import SetupRedis
from ...utils.session import register_model_session
from ..echo import echo


def serve(
    model: str,
    port: int = None,
    track: bool = False,
    tracking_use_case: str = "local",
    enable_local_cache: bool = True,
    local_cache_only: bool = False,
    cloud_cache_only: bool = False,
    cache_only: bool = False,
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

    output_source = None

    if local_cache_only:
        output_source = OutputSource.LOCAL_ONLY
        enable_local_cache = True
    if cloud_cache_only:
        output_source = OutputSource.CLOUD_ONLY
    if cache_only:
        output_source = OutputSource.CACHE_ONLY
        enable_local_cache = True

    mdl = ErsiliaModel(
        model,
        output_source=output_source,
        preferred_port=port,
        cache=enable_local_cache,
        maxmemory=max_cache_memory_frac,
    )

    SetupRedis(enable_local_cache, max_cache_memory_frac)

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
