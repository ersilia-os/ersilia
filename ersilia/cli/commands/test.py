import click

from . import ersilia_cli
from ...publish.test import ModelTester
from ... import ModelBase


def test_cmd():
    """Test a model"""

    # Example usage: ersilia test {MODEL}
    @ersilia_cli.command(
        short_help="Test a model",
        help="Test a model and obtain performance metrics",
    )
    @click.argument("model", type=click.STRING)
    def test(model):
        mdl = ModelBase(model)
        model_id = mdl.model_id
        mt = ModelTester(model_id=model_id)
        click.echo("Checking model information")
        mt.run()
