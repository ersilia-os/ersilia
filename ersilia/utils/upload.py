import os
import subprocess


class OsfUploader(object):
    def __init__(self, overwrite, username, password):
        self.overwrite = overwrite
        self.username = username
        self.password = password

    def push(self, project_id, filename, destination):
        """Upload file to OSF"""
        environ = os.environ.copy()
        environ["OSF_USERNAME"] = self.username
        environ["OSF_PASSWORD"] = self.password
        command = ["osf", "-p", project_id, "upload", filename, destination]
        subprocess.check_call(command, env=environ)
