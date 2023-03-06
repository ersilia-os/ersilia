import click

from . import ersilia_cli
from ...publish.test import LocalModelTester, RemoteModelTester
from ... import ModelBase


def test_cmd():
    """Test a model"""

    # Example usage: ersilia test {MODEL} [--repo_path {DIRECTORY}]
    @ersilia_cli.command(
        short_help="Test a model",
        help="Check that a model will work in the local device",
    )
    @click.argument("model", type=click.STRING)
    @click.option(
        "--repo_path", default=None, help="Local folder where the model is stored"
    )
    def test(model, repo_path):
        mdl = ModelBase(model)
        model_id = mdl.model_id
        if repo_path is not None:
            mt = LocalModelTester(model_id=model_id, repo_path=repo_path)
        else:
            mt = RemoteModelTester(model_id=model_id)
        click.echo("Checking metadata")
        mt.check_metadata()
        click.echo("Checking fetch")
        mt.check_fetch()
