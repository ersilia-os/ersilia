import click
import os

from . import ersilia_cli
from .. import echo
from ...utils.terminal import run_command
from ...utils import tmp_pid_file
from ... import ErsiliaModel


def close_cmd():
    # Example usage: ersilia close {MODEL}
    @ersilia_cli.command(short_help="Close model", help="Close model")
    @click.argument("model", type=click.STRING)
    def close(model):
        mdl = ErsiliaModel(model)
        mdl.close()
        echo(":no_entry: Model {0} closed".format(mdl.model_id), fg="green")
