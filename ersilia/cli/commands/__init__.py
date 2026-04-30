import rich_click as click
from rich_click import RichCommand, RichGroup

from ... import __version__ as __version__
from ... import logger

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


@click.group(cls=ErsiliaCommandGroup, context_settings={"show_default": True}, epilog="To learn more about a specific command, run: ersilia COMMAND --help")
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
