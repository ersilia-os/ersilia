import subprocess
import uuid

class Docker(object):

    def __init__(self):
        self.client = docker.from_env()
        self.tag = "tmp"

    def build(self, path):
        img = str(uuid.uuid4())
        tag = "tmp"
        subprocess.
        self.client.build(path=path, )

    def remove(self, image):
        pass

    def copy(self, )
        pass
