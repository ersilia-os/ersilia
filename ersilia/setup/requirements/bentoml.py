import subprocess
import sys
import tempfile
import os
import shutil


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

        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "version.txt")
        cmd = "bentoml --version > {0}".format(tmp_file)
        subprocess.Popen(cmd, shell=True).wait()
        with open(tmp_file, "r") as f:
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
