import rich_click as click

from ...utils.uninstall import Uninstaller
from .. import echo
from . import ersilia_cli


def uninstall_cmd():
    """Uninstalls all Ersilia artifacts present locally on the user's system"""

    @ersilia_cli.command(
        short_help="Uninstall ersilia",
        help="Fully uninstall Ersilia from the local system. Removes the EOS directory with all fetched models, all Ersilia Docker images, and the ersilia pip package.",
    )
    def uninstall():
        echo(
            "This will remove the EOS directory with all fetched models, all Ersilia Docker images, and the ersilia pip package.",
            fg="yellow",
        )
        if not click.confirm("Are you sure you want to proceed?", default=False):
            echo("Aborted.", fg="green")
            return
        ui = Uninstaller()
        ui.uninstall()
