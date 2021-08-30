import click
from bentoml.cli.click_utils import BentoMLCommandGroup
from ... import __version__
from ... import logger

@click.group(cls=BentoMLCommandGroup)
@click.version_option(version=__version__)
@click.option(
    "--verbose",
    default=False,
    is_flag=True,
    help="Show logging on terminal when running commands."
)
def ersilia_cli(verbose):
    """
    Ersilia CLI
    """
    if verbose:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)
