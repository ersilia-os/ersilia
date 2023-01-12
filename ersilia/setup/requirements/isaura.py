import importlib
from ...utils.terminal import run_command


class IsauraRequirement(object):
    def __init__(self):
        self.name = "isaura"
        try:
            importlib.import_module(self.name)
        except:
            self.install()

    def install(self):
        run_command()
