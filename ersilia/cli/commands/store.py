import click

from ...contrib.store import ModelStorager
from . import ersilia_cli


def store_cmd():
    """Creates store command"""
    # Example usage: ersilia store {MODEL_ID}
    @ersilia_cli.command(
        short_help="Store a model",
        help="Store a model in GitHub and the chosen file storage. "
             "This option is only for developers and requires credentials."
    )
    @click.argument("model_id", type=click.STRING)
    @click.option("--path", required=True, type=click.Path())
    def store(model_id, path):
        ms = ModelStorager()
        ms.store(path, model_id)
