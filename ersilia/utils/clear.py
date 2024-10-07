import os
import shutil

from .conda import SimpleConda
from ..default import EOS, BENTOML_PATH


class Uninstaller(object):
    def __init__(self):
        pass

    def _directories(self):
        shutil.rmtree(EOS)
        shutil.rmtree(BENTOML_PATH)

    def _conda(self):
        sc = SimpleConda()
        for env in sc._env_list():
            if env.startswith("#"):
                continue
            if not env.startswith("eos"):
                continue
            env = env.split(" ")[0]
            if len(env.split("-")[0]) == 7:
                sc.delete(env)

    def uninstall(self):
        self._conda()
        self._directories()

