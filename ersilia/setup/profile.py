from ..default import bashrc_cli_snippet


class SetupProfile(object):
    """
    Class to handle profile setup.
    """

    def __init__(self):
        pass

    def setup(self):
        """
        Set up profile.
        """
        bashrc_cli_snippet()
