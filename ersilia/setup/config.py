import os

from ..default import CONFIG_JSON, EOS
from ..utils.config import Checker


class SetupConfig(object):
    """
    Class to handle configuration setup.
    """

    def __init__(self):
        pass

    def setup(self):
        """
        Set up configuration.
        """
        if self._is_done("config"):
            return
        if os.path.exists(os.path.join(EOS, CONFIG_JSON)):
            return
        checker = Checker()
        checker.config()
