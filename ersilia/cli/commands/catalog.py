import click

from . import ersilia_cli
from ...hub.catalog import ModelCatalog


def catalog_cmd():
    """Creates catalog command"""
    # Example usage: ersilia catalog
    @ersilia_cli.command(
        help="List a catalog of models",
    )
    @click.option(
        '--hub',
        is_flag=True,
        default=False,
        help="Show catalog of models available in the Ersilia Model Hub"
    )
    @click.option(
        '--local',
        is_flag=True,
        default=False,
        help="Show catalog of models available in the local computer"
    )
    @click.option(
       '--backlog',
       is_flag=True,
       default=False,
       help="Show models backlog (wish list) in a Google Spreadsheet"
    )
    def catalog(hub=False, local=False, backlog=False):
        mc = ModelCatalog()
        if hub:
            click.echo(mc.hub())
        elif local:
            click.echo(mc.local())
        elif backlog:
            click.echo(mc.backlog())
        else:
            click.echo(click.style("Specifiy one of --hub --local or --backlog", fg="blue"))
