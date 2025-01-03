import os

import h5py

from ..default import H5_DATA_FILE, ISAURA_GDRIVE, ISAURA_TEAM_GDRIVE
from . import terminal

try:
    from pydrive2.auth import GoogleAuth
    from pydrive2.drive import GoogleDrive
except:
    GoogleAuth = None
    GoogleDrive = None

from .config import Secrets


def set_secrets_file():
    secrets = Secrets()
    if not os.path.exists(secrets.gdrive_client_secrets_json):
        secrets.fetch_gdrive_secrets_from_github()
    GoogleAuth.DEFAULT_SETTINGS["client_config_file"] = (
        secrets.gdrive_client_secrets_json
    )
    return GoogleAuth


class DVCFetcher(object):
    """
    A class to fetch data using DVC (Data Version Control).

    Parameters
    ----------
    local_repo_path : str
        The local repository path.
    """

    def __init__(self, local_repo_path):
        self.repo_path = local_repo_path

    def get_data(self):
        """
        Fetch data using DVC.
        """
        if self.check_dvc_exists():
            terminal.run_command("dvc --cd " + self.repo_path + " pull")

    def check_dvc_exists(self):
        """
        Check if DVC file exists.

        Returns
        -------
        bool
            True if DVC file exists, False otherwise.
        """
        if os.path.isfile(self._data_path() + ".dvc"):
            return True
        return False

    def check_h5_exists(self):
        """
        Check if H5 file exists.

        Returns
        -------
        bool
            True if H5 file exists, False otherwise.
        """
        if os.path.isfile(self._data_path()):
            return True
        return False

    def has_data(self):
        """
        Check if H5 file contains data.

        Returns
        -------
        bool
            True if H5 file contains data, False otherwise.
        """
        if self.check_h5_exists():
            with h5py.File(self._data_path(), "a") as f:
                if len(f.keys()) > 0:
                    return True
        return False

    def _data_path(self):
        return os.path.join(self.repo_path, H5_DATA_FILE)


class DVCBrancher(object):
    """
    A class to manage DVC branches.
    """

    def __init__(self):
        pass


class DVCSetup(object):
    """
    A class to set up DVC with Google Drive.

    Parameters
    ----------
    local_repo_path : str
        The local repository path.
    model_id : str
        The model identifier.
    """

    def __init__(self, local_repo_path, model_id):
        self.repo_path = local_repo_path
        self.model_id = model_id
        GoogleAuth = set_secrets_file()
        gauth = GoogleAuth()
        gauth.LocalWebserverAuth()
        self.drive = GoogleDrive(gauth)

    def gdrive_setup(self):
        """
        Set up Google Drive folder for DVC.
        """
        folder = self.drive.CreateFile(
            {
                "title": self.model_id,
                "parents": [{"id": ISAURA_GDRIVE}],
                "mimeType": "application/vnd.google-apps.folder",
            }
        )
        folder.Upload()

    def gdrive_folder_id(self):
        """
        Get the Google Drive folder ID for the model.

        Returns
        -------
        str
            The Google Drive folder ID.
        """
        fileList = self.drive.ListFile(
            {
                "q": "'" + ISAURA_GDRIVE + "' in parents and trashed=false",
                "corpora": "teamDrive",
                "teamDriveId": ISAURA_TEAM_GDRIVE,
                "includeTeamDriveItems": True,
                "supportsTeamDrives": True,
            }
        ).GetList()

        for file in fileList:
            if file["title"] == self.model_id:
                return str(file["id"])

    def set_dvc_gdrive(self):
        """
        Set up DVC remote storage on Google Drive.
        """
        terminal.run_command("dvc --cd {0} add data.h5".format(self.repo_path))
        terminal.run_command(
            "dvc --cd "
            + self.repo_path
            + " remote add -d public_repo gdrive://"
            + self.gdrive_folder_id()
        )
        cmd = "dvc --cd " + self.repo_path + " push"
        terminal.run_command(cmd, quiet=False)

    def git_add_and_commit(self, message="Set to public data repo"):
        """
        Add and commit changes to Git.

        Parameters
        ----------
        message : str, optional
            The commit message. Default is "Set to public data repo".
        """
        cwd = os.getcwd()
        os.chdir(self.repo_path)
        terminal.run_command(
            "git add data.h5.dvc"
        )  # TODO: Be more specific in the files/folders to be added
        terminal.run_command("git commit -m '{0}'".format(message))
        os.chdir(cwd)
