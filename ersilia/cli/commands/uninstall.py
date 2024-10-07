from . import ersilia_cli
from ...utils.uninstall import Uninstaller


def uninstall_cmd():
    """Uninstalls all contents related to Ersilia available in the local computer"""

    # Example usage: ersilia setup
    @ersilia_cli.command(
        short_help="Uninstalls ersilia",
        help="Uninstalls and removes all contents related to Ersilia available in the local computer.",
    )
    def uninstall():
        cl = Uninstaller()
        cl.uninstall()
