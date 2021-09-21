import click
import os
import json

from . import ersilia_cli
from .. import echo
from .utils.utils import tmp_pid_file
from ...serve.api import Api
from ... import ModelBase


def api_cmd():
    """Create api command"""
    # Example usage: ersilia api {MODEL} {API_NAME} -i {INPUT} [-o {OUTPUT} -b {BATCH_SIZE}]
    @ersilia_cli.command(
        short_help="Run API on a served model", help="Run API on a served model"
    )
    @click.argument("model", type=click.STRING)
    @click.argument("api_name", type=click.STRING)
    @click.option("-i", "--input", "input", required=True, type=click.STRING)
    @click.option(
        "-o", "--output", "output", required=False, default=None, type=click.STRING
    )
    @click.option(
        "-b", "--batch_size", "batch_size", required=False, default=100, type=click.INT
    )
    def api(model, api_name, input, output, batch_size):
        model_id = ModelBase(model).model_id
        tmp_file = tmp_pid_file(model_id)
        if not os.path.exists(tmp_file):
            echo(
                "Model {0} is not served! run 'ersilia serve {0}' first".format(
                    model_id
                ),
                fg="red",
            )
            return
        with open(tmp_file, "r") as f:
            for l in f:
                url = l.rstrip().split()[1]
        api = Api(model_id, url, api_name)
        for result in api.post(input=input, output=output, batch_size=batch_size):
            if result is not None:
                click.echo(json.dumps(result, indent=4))
            else:
                echo(
                    "Something went wrong. Try running in verbose mode (-v) or contact us at hello@ersilia.io",
                    fg="red",
                )
