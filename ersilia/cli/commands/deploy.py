import click

from . import ersilia_cli


def deploy_cmd():
    """Creates deploy command"""
    from ...contrib.deploy import Deployer
    # Example usage: ersilia deploy {MODEL_ID}
    @ersilia_cli.command(
        short_help="Deploy model to the cloud",
        help="Deploy model in a cloud service. "
             "This option is only for developers and requires credentials."
    )
    @click.argument("model_id", type=click.STRING)
    @click.option("--cloud", default="heroku", type=click.STRING)
    def deploy(model_id, cloud):
        dp = Deployer(cloud=cloud)
        if dp.dep is None:
            click.echo(click.style("Please enter a valid cloud option", fg="red"))
            click.echo(click.style("Only 'heroku' and 'local' are available for the moment...", fg="yellow"))
            return
        dp.deploy(model_id)
