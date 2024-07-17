import importlib
import sys
from ...utils.terminal import run_command


class PyAirtableRequirement:
    def __init__(self):
        self.name = "pyairtable"

    def is_installed(self):
        try:
            importlib.import_module(self.name)
            return True
        except:
            return False

    def install(self):
        version = "<2" if sys.version_info.minor == 7 else "<3"
        run_command(f"python -m pip install 'pyairtable{version}'")
