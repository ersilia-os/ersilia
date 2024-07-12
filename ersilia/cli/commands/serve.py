import click
import time

from .. import echo
from . import ersilia_cli
from ... import ErsiliaModel
from ..messages import ModelNotFound
from ...core.tracking import write_persistent_file


def serve_cmd():
    """Creates serve command"""

    # Example usage: ersilia serve {MODEL}
    @ersilia_cli.command(short_help="Serve model", help="Serve model")
    @click.argument("model", type=click.STRING)
    @click.option("--lake/--no-lake", is_flag=True, default=True)
    @click.option("--docker/--no-docker", is_flag=True, default=False)
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
    def serve(model, lake, docker, port, track):
        start_time = time.time()
        if docker:
            service_class = "docker"
        else:
            service_class = None
        mdl = ErsiliaModel(
            model,
            save_to_lake=lake,
            service_class=service_class,
            preferred_port=port,
            track_runs=track,
        )
        if not mdl.is_valid():
            ModelNotFound(mdl).echo()

        mdl.serve()
        if mdl.url is None:
            echo("No URL found. Service unsuccessful.", fg="red")
            return
        echo(
            ":rocket: Serving model {0}: {1}".format(mdl.model_id, mdl.slug), fg="green"
        )
        echo("")
        echo("   URL: {0}".format(mdl.url), fg="yellow")
        echo("   PID: {0}".format(mdl.pid), fg="yellow")
        echo("   SRV: {0}".format(mdl.scl), fg="yellow")
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

        if track:
            """
            Retrieve the time taken in seconds to serve the Model.
            """
            end_time = time.time()
            duration = end_time - start_time
            content = "Total time taken: {0}\n".format(duration)
            write_persistent_file(content, mdl.model_id)
