from . import ersilia_cli


def conda_cmd():
    """Creates conda command"""
    # Example usage: ersilia conda activate eos0aaa
    @ersilia_cli.command(
        short_help="Use conda from ersilia.",
        help="Use conda from ersilia in order to access model environments."
    )
    def conda():
        pass
