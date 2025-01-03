import click

from ...hub.content.card import ModelCard
from ...hub.content.catalog import ModelCatalog
from . import ersilia_cli


def catalog_cmd():
    """
    Creates the catalog command for the CLI.

    This command allows users to list a catalog of models available either locally or in the model hub.
    It provides options to display the catalog in various formats(such as tables by default or json), show more detailed information,
    and view model cards for specific models.

    Returns
    -------
    function
        The catalog command function to be used by the CLI and for testing in the pytest.

    Examples
    --------
    .. code-block:: console

    Display model card for a specific model ID and show catalog in json format:
    $ ersilia catalog --card <model_id> --as-json
    """

    # Example usage: ersilia catalog
    @ersilia_cli.command(help="List a catalog of models")
    @click.option(
        "-h/-l",
        "--hub/--local",
        is_flag=True,
        default=False,
        help="Show catalog of models available in the model hub",
    )
    @click.option(
        "--file_name", "-f", default=None, type=click.STRING, help="Catalog file name"
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
        "-j/-t",
        "--as-json/--as-table",
        is_flag=True,
        default=False,
        help="Show catalog in table format",
    )
    def catalog(
        hub=False,
        file_name=None,
        browser=False,
        more=False,
        card=False,
        model=None,
        as_json=False,
    ):
        if card and not model:
            click.echo(
                click.style("Error: --card option requires a model ID", fg="red"),
                err=True,
            )
            return
        elif card and model:
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
        else:
            mc = ModelCatalog(less=not more)

            if hub:
                catalog_table = mc.hub()
            else:
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
                catalog = (
                    catalog_table.as_json() if as_json else catalog_table.as_table()
                )
            else:
                catalog_table.write(file_name)
                catalog = None

            click.echo(catalog)

    return catalog
