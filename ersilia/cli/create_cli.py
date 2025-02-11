import sys
from ..auth.auth import Auth
from .cmd import Command
from .commands import ersilia_cli

def create_ersilia_cli():
    """
    Creates and configures the Ersilia CLI.

    This function initializes the Command class, checks if the user is a contributor,
    and dynamically imports and executes various CLI commands based on the user's role.

    Returns
    -------
    ersilia_cli : module
        The configured Ersilia CLI module.
    """
    # Check Python version
    if sys.version_info < (3.8,):
        print("\033[91mWARNING: Ersilia does not support Python versions below 3.8. Please upgrade your Python version.\033[0m")

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
    if is_contributor:
        cmd.setup()

    return ersilia_cli
