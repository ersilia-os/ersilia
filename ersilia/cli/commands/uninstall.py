from . import ersilia_cli
from ...utils.uninstall import Uninstaller


def uninstall_cmd():
    """Uninstalls all Ersilia artifacts present locally on the user's system"""

    # Example usage: ersilia setup
    @ersilia_cli.command(
        short_help="Uninstall ersilia",
        help="Uninstalls all Ersilia artifacts present locally on the user's system.",
    )
    def uninstall():
        ui = Uninstaller()
        ui.uninstall()
