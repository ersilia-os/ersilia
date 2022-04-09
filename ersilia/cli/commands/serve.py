import click
import os
import time
import csv

from ...utils import tmp_pid_file
from . import ersilia_cli
from .. import echo
from ... import ErsiliaModel
from ..messages import ModelNotFound


def serve_cmd():
    """Creates serve command"""
    # Example usage: ersilia serve {MODEL}
    @ersilia_cli.command(short_help="Serve model", help="Serve model")
    @click.argument("model", type=click.STRING)
    @click.option("--lake/--no-lake", is_flag=True, default=True)
    @click.option("--docker/--no-docker", is_flag=True, default=False)
    def serve(model, lake, docker):
        if docker:
            service_class = "docker"
        else:
            service_class = None
        mdl = ErsiliaModel(model, save_to_lake=lake, service_class=service_class)
        if not mdl.is_valid():
            ModelNotFound(mdl).echo()
        mdl.serve()
        if mdl.url is None:
            echo("No URL found. Service unsuccessful.", fg="red")
            return
        echo(
            ":rocket: Serving model {0}: {1}".format(mdl.model_id, mdl.slug), fg="green"
        )
        echo("")
        echo("   URL: {0}".format(mdl.url), fg="yellow")
        echo("   PID: {0}".format(mdl.pid), fg="yellow")
        echo("   SRV: {0}".format(mdl.scl), fg="yellow")
        echo("")
        echo(":backhand_index_pointing_right: Available APIs:", fg="blue")
        apis = mdl.get_apis()
        for api in apis:
            echo("   - {0}".format(api), fg="blue")
            
        # update last usage time of model in fetched_models.txt, every time it 'serves' 
        ts_str = str(time.time())
        with open("fetched_models.txt") as infile:
            models = dict(csv.reader(infile))
        infile.close()
        print(mdl.model_id, models.keys())
        if mdl.model_id in models.keys():
            models[mdl.model_id] = ts_str
        with open('fetched_models.txt', 'w') as f:
            for key, values in models.items():
                f.write(f"{key},{values}\n")
