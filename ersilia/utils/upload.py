import os
import subprocess


class OsfUploader(object):
    """
    A class to upload files to the Open Science Framework (OSF).

    Parameters
    ----------
    overwrite : bool
        Whether to overwrite existing files.
    username : str
        The OSF username.
    password : str
        The OSF password.

    Methods
    -------
    push(project_id, filename, destination)
        Upload a file to OSF.
    """

    def __init__(self, overwrite, username, password):
        self.overwrite = overwrite
        self.username = username
        self.password = password

    def push(self, project_id, filename, destination):
        """
        Upload a file to OSF.

        Parameters
        ----------
        project_id : str
            The OSF project ID.
        filename : str
            The path to the file to upload.
        destination : str
            The destination path on OSF.
        """
        environ = os.environ.copy()
        environ["OSF_USERNAME"] = self.username
        environ["OSF_PASSWORD"] = self.password
        command = ["osf", "-p", project_id, "upload", filename, destination]
        subprocess.check_call(command, env=environ)
