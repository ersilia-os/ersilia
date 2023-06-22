import subprocess
import sys


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
        return False

    def install(self):
        print("Installing bentoml (the ersilia version)")
        cmd = "{0} -m pip install -U git+https://github.com/ersilia-os/bentoml-ersilia.git".format(
            sys.executable
        )
        subprocess.Popen(cmd, shell=True).wait()
