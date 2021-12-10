import os

try:
    from isaura.core.hdf5 import Hdf5Explorer
except:
    Hdf5Explorer = None

from ..utils.terminal import run_command

from .. import ErsiliaBase


class IsauraManager(ErsiliaBase):
    def __init__(self, model_id, config_json, credentials_json):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.model_id = model_id
        self.hdf5 = Hdf5Explorer(model_id)

    def remove_local_duplicates(self):
        for api_name in self.hdf5.list_apis():
            self.hdf5.set_curr_api(api_name)
            self.hdf5.remove_local_duplicates()

    def append_local_to_public(self, secret_keys=[]):
        self.pull()
        for api_name in self.hdf5.list_apis():
            self.hdf5.set_curr_api(api_name)
            self.hdf5.append_local_to_public(secret_keys=secret_keys)

    def push(self, message="Add data with DVC"):
        cwd = os.getcwd()
        os.chdir(self._model_path(self.model_id))
        run_command("dvc add data.h5")
        run_command("dvc push")
        run_command("git add data.h5.dvc")
        run_command("git commit -m '{0}'".format(message))
        run_command("git push")
        os.chdir(cwd)

    def pull(self):
        cwd = os.getcwd()
        os.chdir(self._model_path(self.model_id))
        run_command("dvc pull")
        self.remove_local_duplicates()
        os.chdir(cwd)
