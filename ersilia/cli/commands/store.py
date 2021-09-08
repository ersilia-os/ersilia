import click

from ... import ModelBase
from ...contrib.store import ModelStorager
from . import ersilia_cli


def store_cmd():
    """Creates store command"""
    # Example usage: ersilia store {MODEL}
    @ersilia_cli.command(
        short_help="Store a model",
        help="Store a model in GitHub and the chosen file storage. "
        "This option is only for developers and requires credentials.",
    )
    @click.argument("model", type=click.STRING)
    @click.option("--path", required=True, type=click.Path())
    def store(model, path):
        model_id = ModelBase(model).model_id
        ms = ModelStorager()
        ms.store(path, model_id)
