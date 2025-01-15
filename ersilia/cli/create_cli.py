import sys

from ..auth.auth import Auth
from .cmd import Command
from .commands import ersilia_cli
from .messages import VersionNotSupported


def create_ersilia_cli():
    """
    Creates and configures the Ersilia CLI.

    This function initializes the Command class, checks if the user is a contributor,
    and dynamically imports and executes various CLI commands based on the user's role.
    It also verifies that Python 3.8+ is being used.

    Returns
    -------
    ersilia_cli : module | None
        The configured Ersilia CLI module if Python version >= 3.8
        None if the program exits due to version check failure
    Raises
    ------
    SystemExit
        If Python version is less than 3.8
    Notes
    -----
    This function will terminate the program with sys.exit(1) if Python version < 3.8
    """
    # Check Python version
    if sys.version_info.major < 3 or sys.version_info.minor < 8:
        VersionNotSupported(sys.version.major, sys.version.minor)

    is_contributor = Auth().is_contributor()

    cmd = Command()

    cmd.auth()
    cmd.catalog()
    cmd.uninstall()
    cmd.close()
    cmd.delete()
    cmd.example()
    cmd.fetch()
    cmd.info()
    cmd.test()

    # TODO: publishing functionalities
    if is_contributor:
        cmd.publish()

    cmd.serve()
    cmd.run()

    # TODO: functions only for contributors
    # Functions only for contributors
    if is_contributor:
        cmd.setup()

    return ersilia_cli
