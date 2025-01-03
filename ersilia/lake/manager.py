import os

try:
    from isaura.core.hdf5 import Hdf5Explorer
except:
    Hdf5Explorer = None

from .. import ErsiliaBase
from ..utils.terminal import run_command


class IsauraManager(ErsiliaBase):
    """
    Manager class for handling Isaura HDF5 operations.

    Parameters
    ----------
    model_id : str
        Identifier for the model.
    config_json : dict
        Configuration settings in JSON format.
    credentials_json : dict
        Credentials settings in JSON format.

    Attributes
    ----------
    model_id : str
        Identifier for the model.
    hdf5 : Hdf5Explorer
        Instance of Hdf5Explorer for managing HDF5 operations.
    """

    def __init__(self, model_id: str, config_json: dict, credentials_json: dict):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.model_id = model_id
        self.hdf5 = Hdf5Explorer(model_id)

    def remove_local_duplicates(self):
        """
        Remove local duplicate entries in the HDF5 file.
        """
        for api_name in self.hdf5.list_apis():
            self.hdf5.set_curr_api(api_name)
            self.hdf5.remove_local_duplicates()

    def append_local_to_public(self, secret_keys: list = []):
        """
        Append local data to the public repository.

        Parameters
        ----------
        secret_keys : list, optional
            List of secret keys to be used for appending data.
        """
        self.pull()
        for api_name in self.hdf5.list_apis():
            self.hdf5.set_curr_api(api_name)
            self.hdf5.append_local_to_public(secret_keys=secret_keys)

    def push(self, message: str = "Add data with DVC"):
        """
        Push local changes to the remote repository.

        Parameters
        ----------
        message : str, optional
            Commit message for the push operation.
        """
        cwd = os.getcwd()
        os.chdir(self._model_path(self.model_id))
        run_command("dvc add data.h5")
        run_command("dvc push")
        run_command("git add data.h5.dvc")
        run_command("git commit -m '{0}'".format(message))
        run_command("git push")
        os.chdir(cwd)

    def pull(self):
        """
        Pull the latest changes from the remote repository and remove local duplicates.
        """
        cwd = os.getcwd()
        os.chdir(self._model_path(self.model_id))
        run_command("dvc pull")
        self.remove_local_duplicates()
        os.chdir(cwd)
