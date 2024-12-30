import os
import subprocess


class OsfRemover(object):
    """
    A class to remove files from the Open Science Framework (OSF).

    Parameters
    ----------
    username : str
        The OSF username.
    password : str
        The OSF password.
    """

    def __init__(self, username, password):
        self.username = username
        self.password = password

    def remove(self, project_id, filename):
        """
        Remove a file from OSF.

        Parameters
        ----------
        project_id : str
            The OSF project ID.
        filename : str
            The name of the file to remove.
        """
        environ = os.environ.copy()
        environ["OSF_USERNAME"] = self.username
        environ["OSF_PASSWORD"] = self.password
        command = ["osf", "-p", project_id, "remove", filename]
        subprocess.check_call(command, env=environ)
