import sys

import rich_click as click

from ... import ErsiliaModel
from ...core.session import Session
from ...utils.logging import logger
from ...utils.session import register_model_session
from ...utils.terminal import print_serve_summary
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

    @ersilia_cli.command(short_help="Serve model", help="Serve model")
    @click.argument("model", type=click.STRING)
    @click.option(
        "--port",
        "-p",
        default=None,
        type=click.INT,
        help="Preferred port to use (integer)",
    )
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
    @click.option("--access", "-a", default=None)
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
        sess.register_store_status(
            read_store, write_store, access, nearest_neighbors, enable_cache
        )
        store_stat = store_status(read_store, write_store)
        if not is_installed("isaura") and (read_store or write_store):
            echo(
                "Isaura is not installed! Please install isaura in your env by running simply \n>> pip install git+https://github.com/ersilia-os/isaura.git.\nTo start all isaura services, run this command >> isaura engine -s.",
                fg="red",
            )
            logger.error(
                "Isaura is not installed! Please install isaura in your env [pip install git+https://github.com/ersilia-os/isaura.git]! To start all isaura services, execute >> isaura engine -s. "
            )
            sys.exit(1)
        if write_store and access is None:
            echo(
                "You need to specifiy the access as [public or private] to write to store!",
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

        print_serve_summary(
            model_id=mdl.model_id,
            slug=mdl.slug,
            url=mdl.url,
            pid=mdl.pid,
            srv=mdl.scl,
            session_dir=mdl.session._session_dir,
            apis=mdl.get_apis(),
            store_stat=store_stat,
            enable_cache=enable_cache,
            tracking_enabled=bool(track),
            tracking_use_case=tracking_use_case,
        )

        logger.success(f"Model {model} is successfully served!")

    return serve
