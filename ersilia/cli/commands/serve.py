import sys

import click

from ... import ErsiliaModel
from ...store.utils import CacheRetrievingOptions, CacheSavingOptions, OutputSource
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
        default=False,
        help="Enable tracking of this serve session",
    )
    @click.option(
        "--tracking-use-case",
        type=click.Choice(
            ["local", "self-service", "hosted", "test"], case_sensitive=False
        ),
        default="local",
        help="Tracking use case. Options: local, self-service, hosted, test",
    )
    @click.option(
        "--cache-saving",
        "cache_saving",
        type=click.Choice(
            ["disabled", "local", "cloud", "hybrid"], case_sensitive=False
        ),
        default="local",
        help="Cache saving options: disabled, local, cloud, hybrid",
    )
    @click.option(
        "--cache-retrieving",
        "cache_retrieving",
        type=click.Choice(
            ["disabled", "local", "cloud", "hybrid"], case_sensitive=False
        ),
        default=None,
        help="Cache retrieving options: disabled, local, cloud, hybrid",
    )
    @click.option(
        "-co",
        "--cache-only",
        "cache_only",
        is_flag=True,
        default=False,
        help="Only use cache, never fetch new",
    )
    @click.option(
        "--max-cache-memory-frac",
        "max_memory",
        type=click.FLOAT,
        default=None,
        help="Max fraction of memory to dedicate to cache",
    )
    def serve(
        model,
        port,
        track,
        tracking_use_case,
        cache_saving,
        cache_retrieving,
        cache_only,
        max_memory,
    ):
        if cache_saving.lower() in ("cloud", "hybrid"):
            echo(
                f"Warning: cache-saving mode '{cache_saving}' is not supported yet. "
                "Please use 'local' or 'disabled'.",
                fg="yellow",
            )
            sys.exit(1)

        output_source = None
        cache_status = "Disabled"
        enable_local_cache = True

        if (
            cache_retrieving == CacheRetrievingOptions.DISABLED
            and cache_saving == CacheSavingOptions.DISABLED
        ):
            enable_local_cache = False
            cache_only = False

        _status_map = {
            CacheRetrievingOptions.LOCAL: ("Local", OutputSource.LOCAL_ONLY),
            CacheRetrievingOptions.CLOUD: ("Cloud", OutputSource.CLOUD_ONLY),
            CacheRetrievingOptions.HYBRID: ("Hybrid", OutputSource.HYBRID),
        }

        try:
            cache_status, output_source = _status_map[cache_retrieving]
        except KeyError:
            pass

        mdl = ErsiliaModel(
            model,
            output_source=output_source,
            preferred_port=port,
            cache=enable_local_cache,
            maxmemory=max_memory,
            cache_only=cache_only,
            cache_saving_source=cache_saving.lower(),
        )
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
        echo("")
        echo(":person_tipping_hand: Information:", fg="blue")
        echo("   - info", fg="blue")
        echo("")
        echo("🔄 Cache retrieving mode:", fg="blue")
        echo(f"   - {cache_status}", fg="red") if cache_status == "Disabled" else echo(
            f"   - {cache_status}", fg="green"
        )
        echo("")
        echo("💾 Cache Saving mode:", fg="blue")
        echo(
            f"   - {cache_saving.lower().capitalize()}", fg="green"
        ) if cache_saving.lower() != CacheSavingOptions.DISABLED else echo(
            "   - Disabled", fg="red"
        )
        echo("")
        echo(":chart_increasing: Tracking:", fg="blue")
        if track:
            echo("   - Enabled ({0})".format(tracking_use_case), fg="green")
        else:
            echo("   - Disabled", fg="red")

    return serve
