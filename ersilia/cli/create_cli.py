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
    try:
        is_contributor = Auth().is_contributor()
    except Exception:
        # Checking contributor status (a GitHub call) must never break the CLI.
        is_contributor = False

    cmd = Command()

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

    return ersilia_cli
