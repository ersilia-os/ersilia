import click
import csv

from . import ersilia_cli
from .. import echo
from ...hub.delete.delete import ModelFullDeleter
from ... import ModelBase


def delete_cmd():
    """Create delete command"""

    def _delete(md, model_id):
        md.delete(model_id)

    # Example usage: ersilia delete {MODEL}
    @ersilia_cli.command(
        short_help="Delete model from local computer",
        help="Delete model from local computer. The BentoML bundle is deleted, as well as the files stored in "
        "the EOS directory and the Pip-installed package",
    )
    @click.argument("model", type=click.STRING)
    def delete(model):
        model_id = ModelBase(model).model_id
        md = ModelFullDeleter()
        if md.needs_delete(model_id):
            echo("Deleting model {0}".format(model_id))
            _delete(md, model_id)
            echo(
                ":collision: Model {0} deleted successfully!".format(model_id),
                fg="green",
            )
        else:
            echo(
                ":person_tipping_hand: Model {0} is not available locally. No delete is necessary".format(
                    model_id
                ),
                fg="yellow",
            )
