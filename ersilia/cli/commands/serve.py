import click

from ... import ErsiliaModel
from ...store.utils import ModelNotInStore, OutputSource, store_has_model
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
        "--output-source",
        type=click.Choice(OutputSource.ALL),
        default=OutputSource.LOCAL_ONLY,
        required=False,
        help=f"Get outputs from locally hosted model only ({OutputSource.LOCAL_ONLY}), \
            from cloud precalculation store only ({OutputSource.CLOUD_ONLY})",
    )
    @click.option("--lake/--no-lake", is_flag=True, default=True)
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
    @click.option("--cache/--no-cache", is_flag=True, default=True)
    @click.option(
        "--max-cache-memory-frac", "maxmemory", type=click.FLOAT, default=None
    )
    def serve(model, output_source, lake, port, track, cache, maxmemory):
        if OutputSource.is_cloud(output_source):
            if store_has_model(model_id=model):
                echo("Model {0} found in inference store.".format(model))
            else:
                ModelNotInStore(model).echo()
        mdl = ErsiliaModel(
            model,
            output_source=output_source,
            save_to_lake=lake,
            preferred_port=port,
            track_runs=track,
            cache=cache,
            maxmemory=maxmemory,
        )
        redis_setup = SetupRedis(cache, maxmemory)
        if not mdl.is_valid():
            ModelNotFound(mdl).echo()

        mdl.serve()
        if mdl.url is None:
            echo("No URL found. Service unsuccessful.", fg="red")
            return

        register_model_session(mdl.model_id, mdl.session._session_dir)
        echo(
            ":rocket: Serving model {0}: {1}".format(mdl.model_id, mdl.slug), fg="green"
        )
        echo("")
        echo("   URL: {0}".format(mdl.url), fg="yellow")
        echo("   PID: {0}".format(mdl.pid), fg="yellow")
        echo("   SRV: {0}".format(mdl.scl), fg="yellow")
        echo("   Output source: {0}".format(mdl.output_source), fg="yellow")
        echo("")
        echo(":backhand_index_pointing_right: To run model:", fg="blue")
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
        echo(":backhand_index_pointing_right: Caching:", fg="blue")
        echo("   - Enabled", fg="green") if redis_setup._is_amenable()[0] else echo(
            f"   - Disabled: {redis_setup._is_amenable()[1]}", fg="red"
        )

    return serve
