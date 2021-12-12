import click
import os

from . import ersilia_cli
from .. import echo
from ...utils.terminal import run_command
from ...utils import tmp_pid_file
from ... import ErsiliaModel
from ...core.session import Session


def close_cmd():
    # Example usage: ersilia close {MODEL}
    @ersilia_cli.command(short_help="Close model", help="Close model")
    def close():
        model_id = Session(config_json=None).current_model_id()
        if model_id is None:
            echo("No model was served")
            return
        mdl = ErsiliaModel(model_id)
        mdl.close()
        echo(":no_entry: Model {0} closed".format(mdl.model_id), fg="green")
