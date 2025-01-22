import click

from ...utils.installers import base_installer, full_installer
from . import ersilia_cli


def setup_cmd():
    """
    Sets up the environment.

    This command allows users to set up the environment for using the CLI.

    Returns
    -------
    function
        The setup command function to be used by the CLI.

    Examples
    --------
    .. code-block:: console

        Set up the environment with full installation:
        $ ersilia setup --full

        Set up the environment with base installation:
        $ ersilia setup --base
    """

    # Example usage: ersilia setup
    @ersilia_cli.command(
        short_help="Setup ersilia",
        help="Setup ersilia, including building a model-server image, a base environment (eos), rdkit, etc.",
    )
    @click.option(
        "--base",
        is_flag=True,
        default=False,
        help="Install only bare-minimum dependencies.",
    )
    @click.option(
        "--full",
        is_flag=True,
        default=True,
        help="Install all the necessary dependencies.",
    )
    def setup(base=False, full=True):
        if base:
            base_installer()
        elif full:
            full_installer()
        else:
            pass
