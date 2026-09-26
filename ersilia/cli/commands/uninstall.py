from ..echo import confirm, echo
from . import ersilia_cli


def uninstall_cmd():
    """Uninstalls all Ersilia artifacts present locally on the user's system"""

    @ersilia_cli.command(
        short_help="Uninstall Ersilia",
        help="Fully uninstall Ersilia from this computer. Removes the Ersilia folder with all fetched models, all Ersilia Docker images and conda environments, and the ersilia Python package.",
    )
    def uninstall():
        echo(
            "This will remove the Ersilia folder with all fetched models, all Ersilia Docker images and conda environments, and the ersilia Python package.",
            fg="yellow",
        )
        if not confirm("Continue?", default=False):
            echo("Aborted. Nothing was removed.")
            return
        from ...utils.uninstall import Uninstaller

        ui = Uninstaller()
        ui.uninstall()
