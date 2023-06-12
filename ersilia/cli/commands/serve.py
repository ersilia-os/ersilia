import click
from . import ersilia_cli
from .. import echo
from ... import ErsiliaModel
from ..messages import ModelNotFound


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
    def serve(model, lake, docker, port):
        if docker:
            service_class = "docker"
        else:
            service_class = None
        mdl = ErsiliaModel(
            model, save_to_lake=lake, service_class=service_class, preferred_port=port
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
