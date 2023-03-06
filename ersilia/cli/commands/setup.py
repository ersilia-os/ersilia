import click

from . import ersilia_cli
from ...utils.installers import base_installer, full_installer


def setup_cmd():
    """Creates setup command"""

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
