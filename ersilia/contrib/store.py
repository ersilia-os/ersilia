"""
Class used to store model when done developing it.

This functionality is used when developing of a model is done.
"""
import os
import shutil
from .. import ErsiliaBase
from ..utils.zip import Zipper
from ..utils.upload import OsfUploader
from ..utils.remove import OsfRemover


class ModelStorager(ErsiliaBase):
    def __init__(self, config_json=None, credentials_json=None, overwrite=True):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.overwrite = overwrite
        self.osf_up = OsfUploader(
            overwrite=overwrite, username=self.cred.OSF.USER, password=self.cred.OSF.PWD
        )

    def _data_path(self, model_id):
        return os.path.join(self.cred.LOCAL.DATA, model_id)

    def _copy_to_local(self, path, model_id):
        dest_path = self._data_path(model_id)
        if os.path.exists(dest_path):
            if self.overwrite:
                shutil.rmtree(dest_path)
                os.makedirs(dest_path)
            else:
                return
        else:
            os.makedirs(dest_path)
        shutil.copytree(path, os.path.join(dest_path, "model"))

    def _zip(self, model_id):
        data_path = self._data_path(model_id)
        zipper = Zipper(remove=False)
        zip_file = data_path + ".zip"
        zipper.zip(data_path, data_path)
        return zip_file

    def _upload_to_osf(self, model_id):
        zip_file = self._zip(model_id)
        destination = "models/" + model_id + ".zip"
        self.osf_up.push(self.cfg.EXT.OSF_PROJECT, zip_file, destination)
        os.remove(zip_file)

    def store(self, path, model_id):
        """Store model in the local data directory and in OSF."""
        self._copy_to_local(path, model_id)
        self._upload_to_osf(model_id)


class ModelRemover(ErsiliaBase):
    def __init__(self, config_json=None, credentials_json=None):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.osf_rm = OsfRemover(
            username=self.cred.OSF.USER, password=self.cred.OSF.PWD
        )

    def remove(self, model_id):
        zip_file = "models/" + model_id + ".zip"
        self.osf_rm.remove(self.cfg.EXT.OSF_PROJECT, zip_file)
