import click

from . import ersilia_cli
from ...hub.fetch import ModelFetcher


def dockerize_cmd():
    """Create dockerize command"""
    # Example usage: ersilia dockerize {MODEL_ID}
    @ersilia_cli.command(
        short_help="Containerize model using docker",
        help="Containerize model in the BentoML style using docker",
    )
    @click.argument("model_id", type=click.STRING)
    def dockerize(model_id):
        mf = ModelFetcher()
        mf.containerize(model_id=model_id)
