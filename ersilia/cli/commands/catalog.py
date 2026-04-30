import json

import rich_click as click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from ...hub.content.card import ModelCard
from ...hub.content.catalog import ModelCatalog
from . import ersilia_cli

_console = Console()


def _print_catalog(catalog_table):
    from rich.table import Table as RichTable

    col_styles = {
        "Index": ("dim", 6),
        "Identifier": ("bold cyan", 12),
        "Slug": ("green", 28),
        "Title": ("", 36),
        "Task": ("magenta", 16),
        "Output Dimension": ("yellow", 18),
        "Fetched From": ("dim cyan", 14),
    }

    table = RichTable(show_header=True, header_style="bold", border_style="grey50", show_lines=True, expand=False)
    for col in catalog_table.columns:
        style, width = col_styles.get(col, ("", 16))
        table.add_column(col, style=style, min_width=width, no_wrap=col in ("Index", "Identifier", "Slug", "Task"))

    for row in catalog_table.data:
        table.add_row(*[str(v) if v is not None else "" for v in row])

    _console.print(table)


def _print_model_card(metadata_json: str):
    data = json.loads(metadata_json)

    def fmt(value):
        if isinstance(value, list):
            return ", ".join(str(v) for v in value)
        return str(value) if value is not None else "—"

    sections = [
        ("Overview", ["Identifier", "Slug", "Status", "Task", "Subtask"]),
        ("Description", ["Title", "Description", "Interpretation"]),
        ("Input / Output", ["Input", "Input Dimension", "Input Shape", "Output", "Output Dimension", "Output Shape", "Output Type", "Output Consistency"]),
        ("Deployment", ["Deployment", "Source", "Source Type", "Docker Architecture", "DockerHub", "S3"]),
        ("Publication", ["License", "Contributor", "Publication Type", "Publication Year", "Publication", "Source Code"]),
        ("Sizes", ["Model Size", "Environment Size", "Image Size"]),
    ]

    table = Table(show_header=False, box=None, padding=(0, 1), expand=True)
    table.add_column("Field", style="bold cyan", no_wrap=True, min_width=24)
    table.add_column("Value", overflow="fold")

    for section_title, fields in sections:
        table.add_row(Text(f" {section_title}", style="bold magenta on grey15"), "")
        for field in fields:
            if field in data:
                table.add_row(f"  {field}", fmt(data[field]))
        table.add_row("", "")

    title = f"[bold]{data.get('Identifier', '')}[/bold]  ·  {data.get('Title', '')}"
    _console.print(Panel(table, title=title, border_style="cyan"))


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
                if output:
                    if not (output.endswith(".json") or output.endswith(".csv")):
                        click.echo(click.style("Error: output file must have a .json or .csv extension.", fg="red"), err=True)
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
                                writer.writerow([key, value if not isinstance(value, list) else ", ".join(str(v) for v in value)])
                    click.echo(click.style(f"Model card saved to {output}", fg="green"))
                else:
                    _print_model_card(model_metadata)
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
                _print_catalog(catalog_table)
            else:
                if not (output.endswith(".csv") or output.endswith(".json")):
                    click.echo(click.style("Error: output file must have a .csv or .json extension.", fg="red"), err=True)
                    return
                catalog_table.write(output)
                click.echo(click.style(f"Catalog saved to {output}", fg="green"))

    return catalog
