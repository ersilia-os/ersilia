import os
from .identifiers import LongIdentifier
from .terminal import run_command


class SimpleDocker(object):

    def __init__(self):
        self.identifier = LongIdentifier()

    @staticmethod
    def _image_name(org, img, tag):
        return "%s/%s:%s" % (org, img, tag)

    def build(self, path, org, img, tag):
        path = os.path.abspath(path)
        cwd = os.getcwd()
        os.chdir(path)
        cmd = "docker build -t %s %s" % (self._image_name(org, img, tag), path)
        run_command(cmd, quiet=True)
        os.chdir(cwd)

    def remove(self, org, img, tag):
        cmd = "docker rmi %s" % self._image_name(org, img, tag)
        run_command(cmd, quiet=True)

    def run(self, org, img, tag, name):
        if name is None:
            name = self.identifier.encode()
        cmd = "docker run -it -d --name %s %s bash" % (name, self._image_name(org, img, tag))
        run_command(cmd, quiet=True)
        return name

    @staticmethod
    def kill(name):
        cmd = "docker kill %s" % name
        run_command(cmd, quiet=True)

    @staticmethod
    def cp_from_container(name, img_path, local_path):
        local_path = os.path.abspath(local_path)
        cmd = "docker cp %s:%s %s" % (name, img_path, local_path)
        run_command(cmd, quiet=True)

    def cp_from_image(self, img_path, local_path, org, img, tag):
        name = self.run(org, img, tag, name=None)
        self.cp_from_container(name, img_path, local_path)
        self.kill(name)

    @staticmethod
    def exec_container(name, cmd):
        cmd = 'docker exec -i %s bash -c "%s"' % (name, cmd)
        run_command(cmd, quiet=True)

    def exec(self, cmd, org, img, tag, name):
        name = self.run(org, img, tag, name=name)
        self.exec_container(name, cmd)
        self.kill(name)
