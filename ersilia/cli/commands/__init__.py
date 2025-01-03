import functools

import click

from ... import __version__, logger
from ..echo import Silencer


class ErsiliaCommandGroup(click.Group):
    """
    Command group for Ersilia CLI commands.
    """

    NUMBER_OF_COMMON_PARAMS = 2

    @staticmethod
    def bentoml_common_params(func):
        """
        Add common parameters to the command.
        """

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
        """
        Register a new command with common parameters.
        """

        def wrapper(func):
            func = ErsiliaCommandGroup.bentoml_common_params(func)
            func.__click_params__ = (
                func.__click_params__[-self.NUMBER_OF_COMMON_PARAMS :]
                + func.__click_params__[: -self.NUMBER_OF_COMMON_PARAMS]
            )
            return click.Group.command(self, *args, **kwargs)(func)

        return wrapper


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
    """
    🦠 Welcome to Ersilia! 💊
    """
    if verbose:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)
    silencer = Silencer()
    silencer.speak()
    if silent:
        silencer.silence()
