import sys
import subprocess

def verify_setuptools():
    try:
        import setuptools
    except ImportError:
        cmd = f"{sys.executable} -m pip install setuptools"
        result = subprocess.Popen(cmd, shell=True).wait()
        if result!= 0:
            raise RuntimeError("Failed to install setuptools.")
