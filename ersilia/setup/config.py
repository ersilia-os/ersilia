import os

from ..default import EOS, CONFIG_JSON
from ..utils.config import Checker


class SetupConfig(object):
    def __init__(self):
        pass

    def setup(self):
        if self._is_done("config"):
            return
        if os.path.exists(os.path.join(EOS, CONFIG_JSON)):
            return
        checker = Checker()
        checker.config()
