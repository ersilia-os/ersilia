import sys

import rich_click as click

from ..echo import confirm, echo
from . import ersilia_cli


def is_installed(package_name):
    try:
        __import__(package_name)
        return True
    except ImportError:
        return False


def serve_cmd():
    """
    Serves a previously fetched model as a local API.

    The model must already be available locally (see ``ersilia fetch``).
    After serving, use ``ersilia run`` to send predictions to it.

    Returns
    -------
    function
        The serve command function to be used by the CLI and for testing in pytest.

    Examples
    --------
    .. code-block:: console

        Serve a model on a chosen port:
        $ ersilia serve <model_id> --port 8080

        Serve a model with run tracking enabled:
        $ ersilia serve <model_id> --track
    """

    def store_status(read_store, write_store):
        if read_store and not write_store:
            return "Enabled: Read Only"
        if not read_store and write_store:
            return "Enabled: Write Only"
        if read_store and write_store:
            return "Enabled: Read & Write"
        return "Disabled"

    @ersilia_cli.command(
        short_help="Serve a model as a local API",
        help=(
            "Start a local API server for a previously fetched model. "
            "Once served, the model is ready to receive predictions via "
            "`ersilia run`. Optional flags configure session tracking, "
            "in-memory caching, and Isaura store access."
        ),
    )
    @click.argument("model", type=click.STRING)
    @click.option(
        "--port",
        "-p",
        default=None,
        type=click.INT,
        show_default="auto-selected",
        help="Preferred port for the local API server.",
    )
    @click.option(
        "-t",
        "--track",
        "track",
        is_flag=True,
        required=False,
        default=False,
        show_default="off",
        help=(
            "Track model runs (input/output stats, errors, timing, model "
            "metadata) by sending event data to Ersilia's S3 tracking bucket."
        ),
    )
    @click.option(
        "--tracking_use_case",
        type=click.Choice(
            ["local", "self-service", "hosted", "test"], case_sensitive=True
        ),
        required=False,
        default="local",
        hidden=True,
        help=(
            "Deployment context recorded with each tracked run. "
            "One of: local, self-service, hosted, test."
        ),
    )
    @click.option(
        "--enable_cache/--disable_cache",
        is_flag=True,
        default=False,
        show_default=True,
        help=(
            "Enable a Redis in-memory cache for prediction results so repeat "
            "inputs return instantly."
        ),
    )
    @click.option(
        "--read_store",
        "-rs",
        is_flag=True,
        default=False,
        show_default="off",
        help=(
            "Read previously computed predictions from the Isaura store "
            "before running the model. Requires Isaura installed."
        ),
    )
    @click.option(
        "--write_store",
        "-ws",
        is_flag=True,
        default=False,
        show_default="off",
        help=(
            "Write predictions to the Isaura store as they are computed. "
            "Requires Isaura installed and `--access`."
        ),
    )
    @click.option(
        "--access",
        "-a",
        type=click.Choice(["public", "private"]),
        default=None,
        show_default="unset",
        help=(
            "Visibility for predictions written to the Isaura store. "
            "One of: public, private. Required with `--write_store`."
        ),
    )
    @click.option(
        "--nearest_neighbors",
        "--nearest_neigbors",
        "-nn",
        "nearest_neighbors",
        is_flag=True,
        default=False,
        hidden=True,
        help=(
            "When reading from the Isaura store, also return approximate "
            "nearest-neighbor matches for inputs not exactly cached. "
            "This feature is experimental and not yet ready for use."
        ),
    )
    @click.option(
        "--max_cache_memory_frac",
        "max_memory",
        type=click.FloatRange(0, 1, min_open=True),
        default=None,
        show_default="0.3",
        help=("Maximum fraction (0.0-1.0) of system RAM the Redis cache may use."),
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
        from ... import ErsiliaModel
        from ...core.session import Session
        from ...hub.content.information import version_label
        from ...utils.exceptions_utils.exceptions import (
            ModelNotAvailableLocallyError,
        )
        from ...utils.logging import logger
        from ...utils.session import (
            deregister_model_session,
            register_model_session,
        )
        from ...utils.terminal import print_serve_summary
        from ..messages import ModelNotFound

        sess = Session(config_json=None)
        existing_session = sess.get() or {}
        already_served = existing_session.get("model_id")
        if already_served:
            # Check the requested model before touching the one being served,
            # so a typo or an unfetched model does not close it for nothing.
            from ... import ModelBase

            requested = ModelBase(model)
            if not requested.is_available_locally():
                echo(
                    f"Model {requested.model_id} is not available locally.",
                    fg="red",
                )
                echo(f"Fetch it first with 'ersilia fetch {requested.model_id}'.")
                sys.exit(1)
            echo(
                f"Model {already_served} is already being served in this terminal.",
                fg="yellow",
            )
            if not confirm(
                f"Close {already_served} and serve {model} instead?", default=True
            ):
                echo(f"Aborted. {already_served} is still being served.")
                sys.exit(0)
            # Actually close the previously served model so its Docker
            # container and session are shut down before serving the new one.
            previous_service_class = sess.current_service_class()
            previous_mdl = ErsiliaModel(
                already_served, service_class=previous_service_class
            )
            previous_mdl.close()
            deregister_model_session(already_served)
            echo(f"Model {already_served} closed.", fg="green")
        sess.register_store_status(
            read_store, write_store, access, nearest_neighbors, enable_cache
        )
        store_stat = store_status(read_store, write_store)
        if not is_installed("isaura") and (read_store or write_store):
            echo(
                "Isaura is needed to use the store, but it is not installed.", fg="red"
            )
            echo(
                "Install it with 'pip install git+https://github.com/ersilia-os/isaura.git'\nand start its services with 'isaura engine -s'."
            )
            logger.error(
                "Isaura is not installed! Please install isaura in your env [pip install git+https://github.com/ersilia-os/isaura.git]! To start all isaura services, execute >> isaura engine -s. "
            )
            sys.exit(1)
        if write_store and access is None:
            echo("Writing to the store needs an access level.", fg="red")
            echo("Add '--access public' or '--access private'.")
            logger.error("Writing to the store needs --access public or private.")
            sys.exit(1)

        try:
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
        except ModelNotAvailableLocallyError as e:
            echo(f"Model {e.model} is not available locally.", fg="red")
            echo(f"Fetch it first with 'ersilia fetch {e.model}'.")
            sys.exit(1)
        if not mdl.is_valid():
            ModelNotFound(mdl).echo()

        if track:
            track = tracking_use_case
        else:
            track = None

        mdl.serve(track_runs=track)
        if mdl.url is None:
            echo(f"Model {mdl.model_id} could not be started.", fg="red")
            echo("Run the command again with 'ersilia -v serve' to see the details.")
            sys.exit(1)

        register_model_session(mdl.model_id, mdl.session._session_dir)

        info = mdl.info() if mdl.is_valid() else {}
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
            version=version_label(info.get("docker_tag"), info.get("card")),
        )

        echo(
            "Run it with 'ersilia run -i INPUT.csv -o OUTPUT.csv', and stop it with 'ersilia close'."
        )
        logger.success(f"Model {model} is successfully served!")

    return serve
