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
        "--local/--hub",
        is_flag=True,
        default=False,
        help="Show catalog of models available in the local computer",
    )
    @click.option(
        "--file_name", "-f", default=None, type=click.STRING, help="Catalog file name"
    )
    @click.option(
        "--browser", is_flag=True, default=False, help="Show catalog in the browser"
    )
    @click.option(
        "--more/--less",
        is_flag=True,
        default=False,
        help="Show more information than just the EOS identifier",
    )
    def catalog(local=False, file_name=None, browser=False, more=False):
        if local is True and browser is True:
            click.echo(
                "You cannot show the local model catalog in the browser",
                fg="red",
            )
        if more:
            only_identifier = False
        else:
            only_identifier = True
        mc = ModelCatalog(only_identifier=only_identifier)
        if browser:
            mc.airtable()
            return
        catalog_table = mc.local() if local else mc.hub()
        if local and not catalog_table.data:
            click.echo(
                click.style(
                    "No local model is available. Please fetch a model by running 'ersilia fetch' command",
                    fg="red",
                )
            )
            return
        if file_name is None:
            catalog = catalog_table.as_json()
        else:
            catalog_table.write(file_name)
            catalog = None

        click.echo(catalog)
