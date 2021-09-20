import importlib
from ...utils.terminal import run_command


class RdkitRequirement(object):
    def __init__(self):
        self.name = "rdkit"
        try:
            importlib.import_module(self.name)
        except:
            self.install()

    def install(self):
        run_command("pip install rdkit-pypi")
