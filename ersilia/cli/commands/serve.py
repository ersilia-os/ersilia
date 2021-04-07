import click
import os

from ...serve.autoservice import AutoService
from .utils.utils import tmp_pid_file
from . import ersilia_cli
from .. import echo


def serve_cmd():
    """Creates serve command"""
    # Example usage: ersilia serve {MODEL_ID}
    @ersilia_cli.command(
        short_help="Serve model",
        help="Serve model"
    )
    @click.argument("model_id", type=click.STRING)
    def serve(model_id):
        srv = AutoService(model_id)
        srv.serve()
        echo(":rocket: Serving model {0}!".format(model_id), fg="green")
        echo("")
        echo("    URL: {0}".format(srv.service.url), fg="yellow")
        echo("    PID: {0}".format(srv.service.pid), fg="yellow")
        echo("")
        echo("Available APIs:")
        apis=srv.get_apis()
        for api in apis:
            echo("- {0}".format(api))
        tmp_file = tmp_pid_file(model_id)
        with open(tmp_file, "a+") as f:
            f.write("{0} {1}{2}".format(srv.service.pid, srv.service.url, os.linesep))
