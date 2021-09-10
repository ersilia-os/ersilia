import click
import os

from ...serve.autoservice import AutoService
from .utils.utils import tmp_pid_file
from . import ersilia_cli
from .. import echo
from ... import ModelBase


def serve_cmd():
    """Creates serve command"""
    # Example usage: ersilia serve {MODEL}
    @ersilia_cli.command(short_help="Serve model", help="Serve model")
    @click.argument("model", type=click.STRING)
    def serve(model):
        model = ModelBase(model)
        model_id = model.model_id
        slug = model.slug
        srv = AutoService(model_id)
        srv.serve()
        echo(":rocket: Serving model {0}: {1}".format(model_id, slug), fg="green")
        echo("")
        echo("   URL: {0}".format(srv.service.url), fg="yellow")
        echo("   PID: {0}".format(srv.service.pid), fg="yellow")
        echo("")
        echo(":backhand_index_pointing_right: Available APIs:", fg="blue")
        apis = srv.get_apis()
        for api in apis:
            echo("   - {0}".format(api), fg="blue")
        tmp_file = tmp_pid_file(model_id)
        with open(tmp_file, "a+") as f:
            f.write("{0} {1}{2}".format(srv.service.pid, srv.service.url, os.linesep))
