from ...utils.uninstall import Uninstaller
from . import ersilia_cli


def uninstall_cmd():
    """Uninstalls all Ersilia artifacts present locally on the user's system"""
    """
    Uninstalls all Ersilia artifacts present locally on the user's system

    Returns
    -------
    function
        The uninstall command function to be used by the CLI.

    Examples
    --------
    .. code-block:: console

        Uninstall and completely clean
        $ ersilia uninstall --all

        Uninstall and clean only docker images:
        $ ersilia uninstall --docker

        Uninstall and clean only docker images and conda environments:
        $ ersilia uninstall --docker --conda
    """    

    # Example usage: ersilia setup
    @ersilia_cli.command(
        short_help="Uninstall ersilia",
        help="Uninstalls all Ersilia artifacts present locally on the user's system.",
    )
	@click.option(
        "-s",
        "--sessions",
        is_flag=True,
        required=False,
        default=True,
    )
	@click.option(
        "-d",
        "--docker",
        is_flag=True,
        required=False,
        default=True,
    )
	@click.option(
        "-c",
        "--conda",
        is_flag=True,
        required=False,
        default=True,
    )
	@click.option(
        "-a",
        "--all",
        is_flag=True,
        required=False,
        default=True,
    )
    def uninstall(
    	sessions,
    	docker,
    	conda,
    	all
    ):
        ui = Uninstaller()
        ui.uninstall(sessions, docker, conda, all)
