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
        run_command("python -m pip install rdkit-pypi")


class ChemblWebResourceClientRequirement(object):
    def __init__(self):
        self.name = "chembl_webresource_client"
        try:
            importlib.import_module(self.name)
        except:
            self.install()

    def install(self):
        run_command("python -m pip install chembl_webresource_client")
