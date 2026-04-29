import rich_click as click

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

    @ersilia_cli.command(
        short_help="List a catalog of models",
        help="List models available locally or in the Ersilia Model Hub. By default shows locally fetched models in table format. Supports detailed metadata and individual model cards.\n\nFor a full list of models visit [bold cyan][link=https://catalog.ersilia.io/]https://catalog.ersilia.io/[/link][/bold cyan]",
    )
    @click.option(
        "-h/-l",
        "--hub/--local",
        is_flag=True,
        default=False,
        help="--hub lists models available in the Ersilia Model Hub; --local (default) lists locally fetched models.",
    )
    @click.option(
        "--output", "-o", default=None, type=click.STRING, help="Save the catalog to a file. Accepted formats: .csv, .json."
    )
    @click.option(
        "--more/--less",
        is_flag=True,
        default=False,
        help="--more shows additional model metadata; --less (default) shows only the EOS identifier.",
    )
    @click.option(
        "--card",
        is_flag=True,
        default=False,
        help="Display the full model card for a given model ID.",
    )
    @click.option(
        "--task",
        default=None,
        type=click.Choice(["Annotation", "Representation", "Sampling"], case_sensitive=False),
        help="Filter models by task.",
    )
    @click.argument(
        "model",
        type=click.STRING,
        required=False,
    )
    def catalog(
        hub=False,
        output=None,
        browser=False,
        more=False,
        card=False,
        model=None,
        task=None,
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
            mc = ModelCatalog(less=not more, task=task)

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
            if output is None:
                catalog = catalog_table.as_table()
            else:
                if not (output.endswith(".csv") or output.endswith(".json")):
                    click.echo(click.style("Error: output file must have a .csv or .json extension.", fg="red"), err=True)
                    return
                catalog_table.write(output)
                catalog = None

            click.echo(catalog)

    return catalog
