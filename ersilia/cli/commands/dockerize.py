import click

from . import ersilia_cli
from ...hub.fetch.fetch import ModelFetcher
from ... import ModelBase


def dockerize_cmd():
    """Create dockerize command"""
    # Example usage: ersilia dockerize {MODEL}
    @ersilia_cli.command(
        short_help="Containerize model using docker",
        help="Containerize model in the BentoML style using docker",
    )
    @click.argument("model", type=click.STRING)
    def dockerize(model):
        model_id = ModelBase(model).model_id
        mf = ModelFetcher()
        mf.containerize(model_id=model_id)
