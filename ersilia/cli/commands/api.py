import click
import os

from . import ersilia_cli
from .. import echo
from .utils.utils import tmp_pid_file
from ...serve.api import Api
from ... import ModelBase


def api_cmd():
    """Create api command"""
    # Example usage: ersilia api {MODEL} {API_NAME} {INPUT}
    @ersilia_cli.command(
        short_help="Run API on a served model", help="Run API on a served model"
    )
    @click.argument("model", type=click.STRING)
    @click.argument("api_name", type=click.STRING)
    @click.argument("input", type=click.STRING)
    def api(model, api_name, input):
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
        result = api.post(input)
        if result is not None:
            click.echo(result)
        else:
            echo("Something went wrong", fg="red")
