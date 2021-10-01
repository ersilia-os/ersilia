import click
import json

from . import ersilia_cli
from ...hub.content.card import ModelCard
from ...serve.schema import ApiSchema
from ... import ModelBase


def card_cmd():
    """Creates card command"""
    # Example usage: ersilia card {MODEL}
    @ersilia_cli.command(
        short_help="Get model info card",
        help="Get model info card from Ersilia Model Hub.",
    )
    @click.argument("model", type=click.STRING)
    @click.option(
        "-s",
        "--schema",
        is_flag=True,
        default=False,
        help="Show schema of the model",
    )
    def card(model, schema):
        mdl = ModelBase(model)
        model_id = mdl.model_id
        if not schema:
            mc = ModelCard()
            click.echo(mc.get(model_id, as_json=True))
        else:
            ac = ApiSchema(model_id, config_json=None)
            click.echo(json.dumps(ac.get(), indent=4))
