import importlib

from ...utils.terminal import run_command


def verify_setuptools():
    try:
        importlib.import_module("setuptools")
    except ModuleNotFoundError:
        cmd = "python -m pip install setuptools"
        run_command(cmd)
