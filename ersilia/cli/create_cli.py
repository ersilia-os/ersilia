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
    is_contributor = Auth().is_contributor()

    cmd = Command()

    cmd.catalog()
    cmd.uninstall()
    cmd.close()
    cmd.delete()
    cmd.example()
    cmd.fetch()
    cmd.info()
    cmd.test()
    cmd.dump()

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
