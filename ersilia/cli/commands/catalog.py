import click


from . import ersilia_cli
from ...hub.content.catalog import ModelCatalog
from ...hub.content.search import ModelSearcher
from ...hub.content.table_update import table


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
    @click.option(
        "-t",
        "--text",
        default=None,
        type=click.STRING,
        help="Shows the  model related to input keyword",
    )
    @click.option(
        "-m",
        "--mode",
        default=None,
        type=click.STRING,
        help="Shows the  model trained via input mode",
    )
    @click.option(
        "-n", "--next", is_flag=True, default=False, help="Shows the next table"
    )
    @click.option(
        "-p", "--previous", is_flag=True, default=False, help="Shows previous table"
    )
    def catalog(
        local=False, search=None, text=None, mode=None, next=False, previous=False
    ):

        mc = ModelCatalog()
        if not (local or text or mode):
            catalog = mc.hub()
            if not (next or previous):
                catalog = table(catalog).initialise()

            if next:
                catalog = table(catalog).next_table()

            if previous:
                catalog = table(catalog).prev_table()

        if local:
            catalog = mc.local()

        if text:
            catalog = mc.hub()
            catalog = ModelSearcher(catalog).search_text(text)
        if mode:
            catalog = mc.hub()
            catalog = ModelSearcher(catalog).search_mode(mode)
        click.echo(catalog)
