import click

from ... import ErsiliaModel
from ...store.utils import OutputSource
from ...utils.cache import SetupRedis
from ...utils.session import register_model_session
from .. import echo
from ..messages import ModelNotFound
from . import ersilia_cli


def serve_cmd():
    """
    Serves a specified model.

    This command allows users to serve a specified model as an API.

    Returns
    -------
    function
        The serve command function to be used by the CLI and for testing in the pytest.

    Examples
    --------
    .. code-block:: console

        Serve a model by its ID:
        $ ersilia serve <model_id> --port 8080

        Serve a model and track the session:
        $ ersilia serve <model_id> --track
    """

    # Example usage: ersilia serve {MODEL}
    @ersilia_cli.command(short_help="Serve model", help="Serve model")
    @click.argument("model", type=click.STRING)
    @click.option(
        "--port",
        "-p",
        default=None,
        type=click.INT,
        help="Preferred port to use (integer)",
    )
    # Add the new flag for tracking the serve session
    @click.option(
        "-t",
        "--track",
        "track",
        is_flag=True,
        required=False,
        default=False,
    )
    @click.option(
        "--tracking-use-case",
        type=click.Choice(
            ["local", "self-service", "hosted", "test"], case_sensitive=True
        ),
        required=False,
        default="local",
        help="Tracking use case. Options: local, self-service, hosted, test",
    )
    @click.option(
        "--enable-local-cache/--disable-local-cache", is_flag=True, default=True
    )
    @click.option("--local-cache-only", is_flag=True, default=False)
    @click.option("--cloud-cache-only", is_flag=True, default=False)
    @click.option("--cache-only", is_flag=True, default=False)
    @click.option(
        "--max-cache-memory-frac", "max_memory", type=click.FLOAT, default=None
    )
    def serve(
        model,
        port,
        track,
        tracking_use_case,
        enable_local_cache,
        local_cache_only,
        cloud_cache_only,
        cache_only,
        max_memory,
    ):
        output_source = None
        cache_status = "Disabled"
        if local_cache_only:
            output_source = OutputSource.LOCAL_ONLY
            enable_local_cache = True
            cache_status = "Local only"
        if cloud_cache_only:
            output_source = OutputSource.CLOUD_ONLY
            cache_status = "Cloud only"
        if cache_only:
            output_source = OutputSource.CACHE_ONLY
            enable_local_cache = True
            cache_status = "Hybrid (local & cloud)"
        mdl = ErsiliaModel(
            model,
            output_source=output_source,
            preferred_port=port,
            cache=enable_local_cache,
            maxmemory=max_memory,
        )
        redis_setup = SetupRedis(enable_local_cache, max_memory)
        if not mdl.is_valid():
            ModelNotFound(mdl).echo()

        if track:
            track = tracking_use_case
        else:
            track = None

        mdl.serve(track_runs=track)
        if mdl.url is None:
            echo("No URL found. Service unsuccessful.", fg="red")
            return

        register_model_session(mdl.model_id, mdl.session._session_dir)
        echo(
            ":rocket: Serving model {0}: {1}".format(mdl.model_id, mdl.slug), fg="green"
        )
        echo("")
        echo("   URL: {0}".format(mdl.url), fg="yellow")
        if str(mdl.pid) != "-1":
            echo("   PID: {0}".format(mdl.pid), fg="yellow")
        echo("   SRV: {0}".format(mdl.scl), fg="yellow")
        echo("   Session: {0}".format(mdl.session._session_dir), fg="yellow")
        echo("")
        echo(":backhand_index_pointing_right: Run model:", fg="blue")
        echo("   - run", fg="blue")
        apis = mdl.get_apis()
        if apis != ["run"]:
            echo("")
            echo("   These APIs are also valid:", fg="blue")
            for api in apis:
                if api != "run":
                    echo("   - {0}".format(api), fg="blue")
        echo("")
        echo(":person_tipping_hand: Information:", fg="blue")
        echo("   - info", fg="blue")
        echo("")
        echo("🔄 Cache fetching mode:", fg="blue")
        echo(f"   - {cache_status}", fg="red") if cache_status == "Disabled" else echo(
            f"   - {cache_status}", fg="green"
        )
        echo("")
        echo(":floppy_disk: Local cache:", fg="blue")
        echo("   - Enabled", fg="green") if redis_setup._is_amenable()[0] else echo(
            "   - Disabled", fg="red"
        )
        echo("")
        echo(":chart_increasing: Tracking:", fg="blue")
        if track:
            echo("   - Enabled ({0})".format(tracking_use_case), fg="green")
        else:
            echo("   - Disabled", fg="red")

    return serve
