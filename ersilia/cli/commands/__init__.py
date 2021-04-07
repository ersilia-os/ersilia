import click
from bentoml.cli.click_utils import BentoMLCommandGroup
from ... import __version__


@click.group(cls=BentoMLCommandGroup)
@click.version_option(version=__version__)
def ersilia_cli():
    """
    Ersilia CLI tool
    """
