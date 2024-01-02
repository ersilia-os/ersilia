from . import ersilia_cli
from ...utils.clear import Clearer


def clear_cmd():
    """Clears all contents related to Ersilia available in the local computer"""

    # Example usage: ersilia setup
    @ersilia_cli.command(
        short_help="Clear ersilia",
        help="Clears all contents related to Ersilia available in the local computer.",
    )
    def clear():
        cl = Clearer()
        cl.clear()
