import click

from . import ersilia_cli
from ...hub.content.catalog import ModelCatalog
from ...hub.content.search import ModelSearcher


def catalog_cmd():
    """Creates catalog command"""
    # Example usage: ersilia catalog
    @ersilia_cli.command(help="List a catalog of models")
    @click.option(
        "-l",
        "--local",
        is_flag=True,
        default=False,
        help="Show catalog of models available in the local computer",
    )
    @click.option("-s", "--search", default=None, type=click.STRING)
    def catalog(local=False, search=None):
        mc = ModelCatalog()
        if not local:
            catalog = mc.hub()
        else:
            catalog = mc.local()
        if search:
            catalog = ModelSearcher(catalog).search(search)
        click.echo(catalog)
