import click
import json
import types

from . import ersilia_cli
from .. import echo
from ... import ErsiliaModel
from ...core.session import Session


def run_cmd():
    """Create run command"""

    # Example usage: ersilia run -i {INPUT} [-o {OUTPUT} -b {BATCH_SIZE}]
    @ersilia_cli.command(short_help="Run a served model", help="Run a served model")
    @click.option("-i", "--input", "input", required=True, type=click.STRING)
    @click.option(
        "-o", "--output", "output", required=False, default=None, type=click.STRING
    )
    @click.option(
        "-b", "--batch_size", "batch_size", required=False, default=100, type=click.INT
    )
    @click.option(
        "-t/", "--track_run/--no_track_run", "track_run", required=False, default=False
    )
    @click.option(
        "--standard",
        is_flag=True,
        default=False,
        help="Assume that the run is standard and, therefore, do not do so many checks.",
    )
    def run(input, output, batch_size, track_run, standard):
        session = Session(config_json=None)
        model_id = session.current_model_id()
        service_class = session.current_service_class()
        if model_id is None:
            echo(
                "No model seems to be served. Please run 'ersilia serve ...' before.",
                fg="red",
            )
            return
        mdl = ErsiliaModel(
            model_id,
            service_class=service_class,
            config_json=None,
            track_runs=track_run,
        )
        result = mdl.run(
            input=input,
            output=output,
            batch_size=batch_size,
            track_run=track_run,
            try_standard=standard,
        )
        if isinstance(result, types.GeneratorType):
            for result in mdl.run(input=input, output=output, batch_size=batch_size):
                if result is not None:
                    echo(json.dumps(result, indent=4))
                else:
                    echo("Something went wrong", fg="red")
        else:
            echo(result)
