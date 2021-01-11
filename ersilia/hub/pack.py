from ..core.base import ErsiliaBase
from ..utils.docker import SimpleDocker
import os

DEFAULT_TAG = "repo"


class RepoDockerizer(ErsiliaBase):

    def __init__(self, config_json=None, path=None, tag=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        if path is None:
            self.path = os.path.abspath(os.getcwd())
        self.docker = SimpleDocker()
        self.org = self.cfg.HUB.ORG
        if tag is None:
            self.tag = DEFAULT_TAG
        else:
            self.tag = tag

    def _resolve_model_id(self, model_id):
        if model_id is None:
            return os.path.basename(os.path.normpath(self.path))
        return model_id

    @staticmethod
    def _inside_docker():
        with open('/proc/self/cgroup', 'r') as procfile:
            for line in procfile:
                fields = line.strip().split('/')
                if fields[1] == 'docker':
                    return True
        return False

    def _image_name(self, model_id):
        return self.docker._image_name(org=self.org, img=model_id, tag=self.tag)

    def dockerize(self, model_id):
        if self._inside_docker():
            return None
        self.docker.build(path=self.path, org=self.org, img=model_id, tag=self.tag)
        return self._image_name(model_id)


def dockerize(model_id):
    rd = RepoDockerizer()
    return rd.dockerize(model_id)


class RepoCondaCreator(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)


    def create(self, model_id):
        pass


def condacreate(model_id):
    rc = RepoCondaCreator()
    return rcc.create(model_id)
