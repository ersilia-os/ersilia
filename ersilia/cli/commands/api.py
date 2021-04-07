import click
import requests
import os

from . import ersilia_cli
from .. import echo
from .utils.utils import tmp_pid_file


def api_cmd():
    """Create api command"""
    # Example usage: ersilia api {MODEL_ID} {API_NAME} {INPUT}
    @ersilia_cli.command(
        short_help="Run API on a served model",
        help="Run API on a served model"
    )
    @click.argument("model_id", type=click.STRING)
    @click.argument("api_name", type=click.STRING)
    @click.argument("input", type=click.STRING)
    def api(model_id, api_name, input):
        tmp_file = tmp_pid_file(model_id)
        if not os.path.exists(tmp_file):
            echo("Model {0} is not served! run 'ersilia serve {0}' first".format(model_id), fg="red")
            return
        with open(tmp_file, "r") as f:
            for l in f:
                url = l.rstrip().split()[1]
        response = requests.post("{0}/{1}".format(url, api_name), json=input)
        if response.status_code == 200:
            click.echo(response.json())
        else:
            echo("Something went wrong", fg="red")
