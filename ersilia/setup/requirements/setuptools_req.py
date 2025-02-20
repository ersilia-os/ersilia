import importlib
import sys

from ...utils.terminal import run_command


def verify_setuptools():
    try:
        importlib.import_module("setuptools")
    except ModuleNotFoundError:
        cmd = f"{sys.executable} -m pip install setuptools"
        run_command(cmd)
