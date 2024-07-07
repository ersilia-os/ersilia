import click
from . import ersilia_cli
from .. import echo
from ... import ErsiliaModel
from ..messages import ModelNotFound, ModelNotInStore
from ...core.tracking import open_persistent_file
from ...store.api import InferenceStoreApi
from ...store.utils import OutputSource

def serve_cmd():
    """Creates serve command"""

    # Example usage: ersilia serve {MODEL}
    @ersilia_cli.command(short_help="Serve model", help="Serve model")
    @click.argument("model", type=click.STRING)
    @click.option(
        "--output-source",
        type=click.Choice(OutputSource.ALL),
        default=OutputSource.LOCAL_ONLY,
        required=False,
        help=f"Get outputs from locally hosted model only ({OutputSource.LOCAL_ONLY}), \
            from cloud precalculation store only ({OutputSource.CLOUD_ONLY})"
            # or from cloud precalculation store first then locally hosted model for any inputs \
            # that haven't been precalculated ({OutputSource.CLOUD_FIRST})"
    )
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
        "-t/",
        "--track_serve/--no_track_serve",
        "track_serve",
        required=False,
        default=False,
    )
    def serve(model, output_source, lake, docker, port, track_serve):
        if docker:
            service_class = "docker"
        else:
            service_class = None
        if OutputSource.is_cloud(output_source):
            store = InferenceStoreApi(model_id=model)
            if not store.has_model():
                # if output_source == OutputSource.CLOUD_ONLY:
                ModelNotInStore(model).echo()
                # echo(
                #     "Model {0} not found in inference store. Serving model with output-source={1}.".format(model, OutputSource.LOCAL_ONLY),
                #     fg="yellow"
                # )
                # output_source = OutputSource.LOCAL_ONLY
            else:
                echo("Model {0} found in inference store.".format(model))
        mdl = ErsiliaModel(
            model, output_source=output_source, save_to_lake=lake, service_class=service_class, preferred_port=port
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

        # Setup persistent tracking
        if track_serve:
            open_persistent_file(mdl.model_id)
