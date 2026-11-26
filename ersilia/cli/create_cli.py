from .cmd import Command
from .commands import ersilia_cli


def create_ersilia_cli():
    """
    Creates and configures the Ersilia CLI.

    This function initializes the Command class and registers the available
    CLI commands.

    Returns
    -------
    ersilia_cli : module
        The configured Ersilia CLI module.
    """
    cmd = Command()

    cmd.catalog()
    cmd.uninstall()
    cmd.close()
    cmd.delete()
    cmd.example()
    cmd.fetch()
    cmd.info()
    cmd.test()
    cmd.serve()
    cmd.run()

    return ersilia_cli
