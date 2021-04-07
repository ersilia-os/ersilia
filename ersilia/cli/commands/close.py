import click
import os

from . import ersilia_cli
from .. import echo
from ...utils.terminal import run_command
from .utils.utils import tmp_pid_file


def close_cmd():
    # Example usage: ersilia close {MODEL_ID}
    @ersilia_cli.command(
        short_help="Close model",
        help="Close model"
    )
    @click.argument("model_id", type=click.STRING)
    def close(model_id):
        tmp_file = tmp_pid_file(model_id)
        with open(tmp_file, "r") as f:
            for l in f:
                pid = int(l.rstrip().split()[0])
                cmd = "kill {0}".format(pid)
                run_command(cmd, quiet=True)
        os.remove(tmp_file)
        echo(":no_entry: Model {0} closed".format(model_id), fg="green")
