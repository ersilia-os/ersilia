from . import terminal
import h5py
import os
from ..default import H5_DATA_FILE, ISAURA_GDRIVE
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive


class DVCFetcher(object):
    def __init__(self, local_repo_path):
        self.repo_path = local_repo_path

    def get_data(self):
        terminal.run_command("dvc --cd " + self.repo_path + " pull")

    def check_dvc_exists(self):
        if os.path.isfile(self._data_path() + ".dvc"):
            return True
        return False

    def check_h5_exists(self):
        if os.path.isfile(self._data_path()):
            return True
        return False

    def has_data(self):
        if self.check_h5_exists():
            with h5py.File(self._data_path(), "a") as f:
                if len(f.keys()) > 0:
                    return True
        return False

    def _data_path(self):
        return os.path.join(self.repo_path, H5_DATA_FILE)


class DVCBrancher(object):
    def __init__(self):
        pass


class DVCSetup(object):
    def __init__(self, local_repo_path, model_id):
        self.repo_path = local_repo_path
        self.model_id = model_id
        gauth = GoogleAuth().LocalWebserverAuth()
        drive = GoogleDrive(gauth)

    def gdrive_setup(self):
        folder = drive.CreateFile(
            {
                "title": self.model_id,
                "parents": [{"id": ISAURA_GDRIVE}],
                "mimeType": "application/vnd.google-apps.folder"
            }
        )
        folder.Upload()

    def gdrive_folder_id(self):
        folder = self.drive.ListFile(
            {
                "q": "title='"
                + self.model_id
                + "' and mimeType='application/vnd.google-apps.folder' and trashed=false"
            }
        ).GetList()
        return folder["id"]

    def set_dvc_gdrive(self):
        terminal.run_command(
            "dvc --cd "
            + self.repo_path
            + " remote add -d public_repo gdrive://"
            + self.gdrive_folder_id()
        )
        terminal.run_command("dvc --cd " + self.repo_path + " push")
