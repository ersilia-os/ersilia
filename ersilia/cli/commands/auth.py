from . import ersilia_cli


def auth_cmd():
    """Creates auth command"""

    # Example usage: ersilia auth login
    @ersilia_cli.command(
        short_help="Log in to ersilia to enter contributor mode.",
        help="Log in to ersilia to enter contributor mode. "
        "GitHub credentials are used.",
    )
    def auth():
        """
        In the user's system profile there is a redirect to gh auth (GitHub authorisation CLI)
        """
        pass
