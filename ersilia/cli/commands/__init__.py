from pathlib import Path

import rich_click as click
from rich_click import RichCommand, RichGroup

from ... import __version__ as __version__
from ... import logger

_ASCII_ART_FILE = Path(__file__).parents[3] / "assets" / "ascii-art.txt"


def _print_logo():
    if not _ASCII_ART_FILE.exists():
        return
    from rich.console import Console
    from rich.text import Text

    # Ersilia palette — deeper tones visible on both dark and light terminals
    gradient = ["#E07040", "#C060C0", "#8855BB", "#5599DD", "#55BB77"]

    console = Console(highlight=False)
    lines = _ASCII_ART_FILE.read_text().splitlines()
    non_empty = [l for l in lines if l.strip()]
    n = max(len(non_empty) - 1, 1)
    ci = 0
    for line in lines:
        if line.strip():
            color = gradient[ci * (len(gradient) - 1) // n]
            console.print(Text(line, style=f"bold {color}"), end="\n")
            ci += 1
        else:
            console.print("")

click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.COLOR_SYSTEM = "truecolor"
click.rich_click.STYLE_OPTION = "bold magenta"
click.rich_click.STYLE_COMMAND = "bold green"
click.rich_click.STYLE_METAVAR = "italic yellow"
click.rich_click.STYLE_SWITCH = "underline cyan"
click.rich_click.STYLE_USAGE = "bold blue"
click.rich_click.STYLE_OPTION_DEFAULT = "dim italic"
click.rich_click.STYLE_HELPTEXT = ""
click.rich_click.STYLE_ERRORS_SUGGESTION = "bold"
click.rich_click.HEADER_TEXT = f"Ersilia version: {__version__}"

# ruff: noqa: D101, D102


class ErsiliaCommandGroup(RichGroup):
    def command(self, *args, **kwargs):
        kwargs.setdefault("cls", RichCommand)
        return RichGroup.command(self, *args, **kwargs)

    def main(self, *args, **kwargs):
        import sys
        argv = sys.argv[1:]
        if not argv or "--help" in argv or "-h" in argv:
            _print_logo()
        return super().main(*args, **kwargs)


@click.group(
    cls=ErsiliaCommandGroup,
    context_settings={"show_default": True},
    epilog="To learn more about a specific command, run: ersilia COMMAND --help",
)
@click.version_option(version=__version__)
@click.option(
    "-v",
    "--verbose",
    default=False,
    is_flag=True,
    help="Show logging on terminal when running commands.",
)
def ersilia_cli(verbose):
    if verbose:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)
