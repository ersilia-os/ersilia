import json

import rich_click as click

from .. import echo
from . import ersilia_cli


def _print_catalog(catalog_table):
    from rich import box
    from rich.padding import Padding
    from rich.table import Table as RichTable

    from ...utils.echo import console

    # Same look as the other panels: default colour, secondary columns dimmed,
    # a quiet border, indented to line up with the icons of other lines.
    secondary = {"Index", "Task", "Output Dimension", "Fetched From"}
    table = RichTable(
        box=box.SIMPLE_HEAD,
        header_style="dim",
        border_style="bright_black",
        expand=False,
        show_edge=False,
        pad_edge=False,
    )
    for col in catalog_table.columns:
        table.add_column(
            col,
            style="dim" if col in secondary else None,
            no_wrap=col in ("Index", "Identifier", "Slug", "Task"),
        )
    for row in catalog_table.data:
        table.add_row(*[str(v) if v is not None else "" for v in row])
    console.print(Padding(table, (0, 0, 0, 2), expand=False))


def _print_model_card(metadata_json: str):
    from ...hub.content.information import print_card_panel

    print_card_panel(json.loads(metadata_json))


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
        "--output",
        "-o",
        default=None,
        type=click.STRING,
        help="Save the catalog to a file. Accepted formats: .csv, .json.",
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
        type=click.Choice(
            ["Annotation", "Representation", "Sampling"], case_sensitive=False
        ),
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
        from ...hub.content.card import ModelCard
        from ...hub.content.catalog import ModelCatalog
        from ..messages import wrong_extension

        if card and not model:
            echo(
                "The --card option needs a model, e.g. 'ersilia catalog --card eos42ez'.",
                fg="red",
                err=True,
            )
            return
        elif card and model:
            try:
                mc = ModelCard()
                model_metadata = mc.get(model, as_json=True)

                if not model_metadata:
                    echo(
                        f"No information was found for model {model}.",
                        fg="red",
                        err=True,
                    )
                    return
                if output:
                    if not output.endswith((".json", ".csv")):
                        wrong_extension([".json", ".csv"], err=True)
                        return
                    data = json.loads(model_metadata)
                    if output.endswith(".json"):
                        with open(output, "w") as f:
                            f.write(model_metadata)
                    else:
                        import csv

                        with open(output, "w", newline="") as f:
                            writer = csv.writer(f)
                            writer.writerow(["Field", "Value"])
                            for key, value in data.items():
                                writer.writerow(
                                    [
                                        key,
                                        value
                                        if not isinstance(value, list)
                                        else ", ".join(str(v) for v in value),
                                    ]
                                )
                    echo(f"Model card saved to {output}.", fg="green")
                else:
                    _print_model_card(model_metadata)
            except Exception as e:
                echo(
                    f"Could not get the information of model {model}: {e}",
                    fg="red",
                    err=True,
                )
            return
        else:
            mc = ModelCatalog(less=not more, task=task)

            if hub:
                catalog_table = mc.hub()
            else:
                catalog_table = mc.local()
                if not catalog_table.data:
                    echo("No models are available locally.", fg="yellow")
                    echo("Fetch one with 'ersilia fetch MODEL'.")
                    return
            if output is None:
                _print_catalog(catalog_table)
            else:
                if not output.endswith((".json", ".csv")):
                    wrong_extension([".json", ".csv"], err=True)
                    return
                catalog_table.write(output)
                echo(f"Catalog saved to {output}.", fg="green")

    return catalog
