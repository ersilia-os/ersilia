import click

from . import ersilia_cli
from ...contrib.deploy import Deployer
from ... import ModelBase


def deploy_cmd():
    """Creates deploy command"""

    # Example usage: ersilia deploy {MODEL}
    @ersilia_cli.command(
        short_help="Deploy model to the cloud",
        help="Deploy model in a cloud service. "
        "This option is only for developers and requires credentials.",
    )
    @click.argument("model", type=click.STRING)
    @click.option("--cloud", default="heroku", type=click.STRING)
    def deploy(model, cloud):
        model_id = ModelBase(model).model_id
        dp = Deployer(cloud=cloud)
        if dp.dep is None:
            click.echo(click.style("Please enter a valid cloud option", fg="red"))
            click.echo(
                click.style(
                    "Only 'heroku' and 'local' are available for the moment...",
                    fg="yellow",
                )
            )
            return
        dp.deploy(model_id)
