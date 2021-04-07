import click

from . import ersilia_cli
from ...hub.card import ModelCard


def card_cmd():
    """Creates card command"""
    # Example usage: ersilia card {MODEL_ID}
    @ersilia_cli.command(
        short_help="Get model info card",
        help="Get model info card from Ersilia Model Hub."
    )
    @click.argument("model_id", type=click.STRING)
    def card(model_id):
        mc = ModelCard()
        click.echo(mc.get(model_id, as_json=True))
