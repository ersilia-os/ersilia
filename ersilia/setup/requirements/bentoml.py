import subprocess
import sys
import os

from ...default import EOS

class BentoMLRequirement(object):
    def __init__(self):
        pass

    def is_installed(self):
        try:
            import bentoml

            return True
        except ImportError:
            return False

    def is_bentoml_ersilia_version(self):
        if not self.is_installed():
            return False

        version_file = os.path.join(EOS, "bentomlversion.txt")
        
        if not os.path.exists(version_file):
            cmd = "bentoml --version > {0}".format(version_file)
            subprocess.Popen(cmd, shell=True).wait()

        with open(version_file, "r") as f:
            text = f.read()
        if "0.11.0" in text:
            return True
        else:
            return False

    def install(self):
        print("Installing bentoml (the ersilia version)")
        cmd = "{0} -m pip install -U git+https://github.com/ersilia-os/bentoml-ersilia.git".format(
            sys.executable
        )
        subprocess.Popen(cmd, shell=True).wait()
