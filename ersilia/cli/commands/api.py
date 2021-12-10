import click
import os
import json
import types

from . import ersilia_cli
from .. import echo
from ... import ErsiliaModel


def api_cmd():
    """Create api command"""
    # Example usage: ersilia api {MODEL} {API_NAME} -i {INPUT} [-o {OUTPUT} -b {BATCH_SIZE}]
    @ersilia_cli.command(
        short_help="Run API on a served model", help="Run API on a served model"
    )
    @click.argument("model", default=None, type=click.STRING)
    @click.argument("api_name", required=False, default=None, type=click.STRING)
    @click.option("-i", "--input", "input", required=True, type=click.STRING)
    @click.option(
        "-o", "--output", "output", required=False, default=None, type=click.STRING
    )
    @click.option(
        "-b", "--batch_size", "batch_size", required=False, default=100, type=click.INT
    )
    def api(model, api_name, input, output, batch_size):
        mdl = ErsiliaModel(model)
        result = mdl.api(
            api_name=api_name, input=input, output=output, batch_size=batch_size
        )
        if isinstance(result, types.GeneratorType):
            for result in mdl.api(
                api_name=api_name, input=input, output=output, batch_size=batch_size
            ):
                if result is not None:
                    click.echo(json.dumps(result, indent=4))
                else:
                    click.echo("Something went wrong", fg="red")
        else:
            click.echo(result)
