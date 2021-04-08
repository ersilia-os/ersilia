import os
import subprocess


class OsfRemover(object):
    def __init__(self, username, password):
        self.username = username
        self.password = password

    def remove(self, project_id, filename):
        """Remove file from OSF"""
        environ = os.environ.copy()
        environ["OSF_USERNAME"] = self.username
        environ["OSF_PASSWORD"] = self.password
        command = ["osf", "-p", project_id, "remove", filename]
        subprocess.check_call(command, env=environ)
