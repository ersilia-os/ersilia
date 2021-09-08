import click

from . import ersilia_cli
from ...hub.content.card import ModelCard
from ... import ModelBase


def card_cmd():
    """Creates card command"""
    # Example usage: ersilia card {MODEL}
    @ersilia_cli.command(
        short_help="Get model info card",
        help="Get model info card from Ersilia Model Hub.",
    )
    @click.argument("model", type=click.STRING)
    def card(model):
        mdl = ModelBase(model)
        model_id = mdl.model_id
        mc = ModelCard()
        click.echo(mc.get(model_id, as_json=True))
