import click

from . import ersilia_cli
from ...hub.content.card import ModelCard, LakeCard
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
        "-s", "--schema", is_flag=True, default=False, help="Show schema of the model"
    )
    @click.option(
        "-l",
        "--lake",
        is_flag=True,
        default=False,
        help="Show the properties of the data lake",
    )
    @click.option(
        "-a", "--api", is_flag=True, default=False, help="Show a list of available API"
    )
    def card(model, schema, lake, api):
        mdl = ModelBase(model)
        model_id = mdl.model_id
        if schema:
            mc = ModelCard()
            click.echo(mc.get(model_id, as_json=True))  # TODO
            return
        if lake:
            mc = LakeCard()
            click.echo(mc.get(model_id, as_json=True))
            return
        if api:
            pass  # TODO
        mc = ModelCard()
        click.echo(mc.get(model_id, as_json=True))
