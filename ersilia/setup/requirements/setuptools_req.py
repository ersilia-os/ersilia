import importlib
import subprocess
import sys


def verify_setuptools():
    try:
        importlib.import_module("setuptools")
    except ModuleNotFoundError:
        cmd = f"{sys.executable} -m pip install setuptools"
        result = subprocess.Popen(cmd, shell=True).wait()
        if result != 0:
            raise RuntimeError("Failed to install setuptools.")
