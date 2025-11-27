import functools

import rich_click as click
from rich_click import RichCommand, RichGroup

from ... import __version__ as __version__
from ... import logger
from ..echo import Silencer

click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.COLOR_SYSTEM = "truecolor"
click.rich_click.STYLE_OPTION = "bold magenta"
click.rich_click.STYLE_COMMAND = "bold green"
click.rich_click.STYLE_METAVAR = "italic yellow"
click.rich_click.STYLE_SWITCH = "underline cyan"
click.rich_click.STYLE_USAGE = "bold blue"
click.rich_click.STYLE_OPTION_DEFAULT = "dim italic"

# ruff: noqa: D101, D102


class ErsiliaCommandGroup(RichGroup):
    NUMBER_OF_COMMON_PARAMS = 2

    @staticmethod
    def bentoml_common_params(func):
        @click.option(
            "-q",
            "--quiet",
            is_flag=True,
            default=False,
            help="Hide all warnings and info logs",
        )
        @functools.wraps(func)
        def wrapper(quiet, *args, **kwargs):
            return func(*args, **kwargs)

        return wrapper

    def command(self, *args, **kwargs):
        kwargs.setdefault("cls", RichCommand)

        def decorator(func):
            func = ErsiliaCommandGroup.bentoml_common_params(func)
            func.__click_params__ = (
                func.__click_params__[-self.NUMBER_OF_COMMON_PARAMS :]
                + func.__click_params__[: -self.NUMBER_OF_COMMON_PARAMS]
            )
            return RichGroup.command(self, *args, **kwargs)(func)

        return decorator


@click.group(cls=ErsiliaCommandGroup)
@click.version_option(version=__version__)
@click.option(
    "-v",
    "--verbose",
    default=False,
    is_flag=True,
    help="Show logging on terminal when running commands.",
)
@click.option(
    "-s",
    "--silent",
    default=False,
    is_flag=True,
    help="Do not echo any progress message.",
)
def ersilia_cli(verbose, silent):
    if verbose:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)

    silencer = Silencer()
    silencer.speak()
    if silent:
        silencer.silence()
