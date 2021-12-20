import tempfile
import os
import shutil

from .. import ErsiliaBase

from .terminal import run_command


class SimpleVenv(ErsiliaBase):
    def __init__(self, root):
        ErsiliaBase.__init__(self, config_json=None, credentials_json=None)
        self.root = os.path.abspath(root)
        self.logger.debug("Setting virtual environment at {0}".format(self.root))

    def _get_path(self, environment):
        return os.path.join(self.root, environment)

    def exists(self, environment):
        if os.path.exists(self._get_path(environment)):
            return True
        else:
            return False

    def create(self, environment):
        path = self._get_path(environment)
        if self.exists(path):
            return
        run_command("python -m venv {0} --symlinks --clear".format(path))

    def delete(self, environment):
        path = self._get_path(environment)
        if not self.exists(path):
            return
        shutil.rmtree(path)

    def run_commandlines(self, environment, commandlines):
        if not self.exists(environment):
            raise Exception("{0} environment does not exist".format(environment))
        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        tmp_script = os.path.join(tmp_folder, "script.sh")
        with open(tmp_script, "w") as f:
            f.write("cd {0}{1}".format(self.root, os.linesep))
            f.write("source {0}/bin/activate{1}".format(environment, os.linesep))
            f.write("{0}{1}".format(commandlines, os.linesep))
            f.write("deactivate")
        run_command("bash {0}".format(tmp_script))
