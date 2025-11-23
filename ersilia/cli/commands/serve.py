import sys

import rich_click as click

from ... import ErsiliaModel
from ...core.session import Session
from ...utils.logging import logger
from ...utils.session import register_model_session
from .. import echo
from ..messages import ModelNotFound
from . import ersilia_cli


def is_installed(package_name):
    try:
        __import__(package_name)
        return True
    except ImportError:
        return False


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

    def store_status(read_store, write_store):
        if read_store and not write_store:
            return "Enabled: Read Only"
        if not read_store and write_store:
            return "Enabled: Writter Only"
        if read_store and write_store:
            return "Enabled: Reader & Writter "
        return "Disabled"

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
    @click.option("--enable-cache/--disable-cache", is_flag=True, default=False)
    @click.option("--read-store", "-rs", is_flag=True, default=False)
    @click.option("--write-store", "-ws", is_flag=True, default=False)
    @click.option("--access", default=None)
    @click.option(
        "--nearest-neigbors", "-nn", "nearest_neighbors", is_flag=True, default=False
    )
    @click.option(
        "--max-cache-memory-frac", "max_memory", type=click.FLOAT, default=None
    )
    def serve(
        model,
        port,
        track,
        tracking_use_case,
        enable_cache,
        read_store,
        write_store,
        access,
        nearest_neighbors,
        max_memory,
    ):
        sess = Session(config_json=None)
        sess.register_store_status(read_store, write_store, access, nearest_neighbors)
        store_stat = store_status(read_store, write_store)
        if not is_installed("isaura") and read_store:
            if not logger.verbosity:
                echo(
                    "Isaura is not installed! Please install isaura in your env by running simply \n>> pip install git+https://github.com/ersilia-os/isaura.git.\nTo start all isaura services, run this command >> isaura engine -s.",
                    fg="red",
                )
            logger.error(
                "Isaura is not installed! Please install isaura in your env [pip install git+https://github.com/ersilia-os/isaura.git]! To start all isaura services, execute >> isaura engine -s. "
            )
            sys.exit(1)
        if write_store and access is None or read_store and access is None:
            if not logger.verbosity:
                echo(
                    "You need to specifiy the access as [public or private] to read/write to store!",
                    fg="red",
                )
            logger.error(
                "You need to specifiy the access as [public or private] to write to store!"
            )
            sys.exit(1)

        mdl = ErsiliaModel(
            model,
            output_source=None,
            preferred_port=port,
            cache=enable_cache,
            maxmemory=max_memory,
            read_store=read_store,
            write_store=write_store,
            access=access,
            nearest_neighbors=nearest_neighbors,
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
        echo("🏬 Isaura Store:", fg="blue")
        echo(f"   - {store_stat}", fg="red") if store_stat == "Disabled" else echo(
            f"   - {store_stat}", fg="green"
        )
        echo("")
        echo(":floppy_disk: Local cache:", fg="blue")
        echo("   - Enabled", fg="green") if enable_cache else echo(
            "   - Disabled", fg="red"
        )
        echo("")
        echo(":chart_increasing: Tracking:", fg="blue")
        if track:
            echo("   - Enabled ({0})".format(tracking_use_case), fg="green")
        else:
            echo("   - Disabled", fg="red")
        logger.success(f"Model {model} is successfully served!")

    return serve
