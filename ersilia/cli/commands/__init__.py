import click
from bentoml.cli.click_utils import BentoMLCommandGroup
from ... import __version__
from ... import logger
from ..echo import Silencer


@click.group(cls=BentoMLCommandGroup)
@click.version_option(version=__version__)
@click.option(
    "-v",
    "--verbose",
    default=False,
    is_flag=True,
    help="Show logging on terminal when running commands.",
)
@click.option(
    "-s",
    "--silent",
    default=False,
    is_flag=True,
    help="Do not echo any progress message.",
)
def ersilia_cli(verbose, silent):
    """
    Ersilia CLI
    """
    if verbose:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)
    silencer = Silencer()
    silencer.speak()  # To reset default
    if silent:
        silencer.silence()
