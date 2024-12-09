import click

from . import ersilia_cli
from ...hub.content.catalog import ModelCatalog
from ...hub.content.card import ModelCard


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
        "-h",
        "--hub",
        is_flag=True,
        default=False,
        help="Show catalog of models available in the model hub",
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
    @click.option(
        "--card",
        is_flag=True,
        default=False,
        help="Use this flag to display model card for a given model ID",
    )
    @click.argument(
        "model",
        type=click.STRING,
        required=False,
    )
    @click.option(
        "--as-table",
        is_flag=True,
        default=False,
        help="Show catalog in table format",
    )
    def catalog(
        local=False,
        hub=False,
        file_name=None,
        browser=False,
        more=False,
        card=False,
        model=None,
        as_table=False,
    ):
        if card and not model:
            click.echo(
                click.style("Error: --card option requires a model ID", fg="red"),
                err=True,
            )
            return
        if card and model:
            try:
                mc = ModelCard()
                model_metadata = mc.get(model, as_json=True)

                if not model_metadata:
                    click.echo(
                        click.style(
                            f"Error: No metadata found for model ID '{model}'", fg="red"
                        ),
                        err=True,
                    )
                    return
                click.echo(model_metadata)
            except Exception as e:
                click.echo(click.style(f"Error fetching model metadata: {e}", fg="red"))
            return
        
        # The idea here is to deter the user from running ersilia catalog --local --hub
        if local and hub:
            click.echo(
                click.style(
                    "Error: Cannot show local and hub models together", fg="red"
                ),
                err=True,
            )
            return
        
        mc = ModelCatalog()
        mc.only_identifier = False if more else True

        if hub:
            if browser:
                mc.airtable()
                return
            
            catalog_table = mc.hub()

        else: # This will work even if the user doesn't explicitly specify the --local flag
            if browser:
                click.echo(
                    click.style(
                        "Error: Cannot show local models in the browser.\nPlease use the --hub option to see models in the browser.",
                        fg="red"
                    )
                )
                return
            catalog_table = mc.local()
            if not catalog_table.data:
                click.echo(
                click.style(
                    "No local models available. Please fetch a model by running 'ersilia fetch' command",
                    fg="red",
                    )
                )
                return
        
        if file_name is None:
            catalog = catalog_table.as_table() if as_table else catalog_table.as_json()
        else:
            catalog_table.write(file_name)
            catalog = None

        click.echo(catalog)

    return catalog
