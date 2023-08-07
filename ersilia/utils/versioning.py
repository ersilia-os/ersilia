import sys
import os

try:
    from bentoml import __version__ as __bentoml_version__
except:
    __bentoml_version__ = None
from .. import ErsiliaBase
from .. import __version__ as __ersilia_version__


class Versioner(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def python_version(self, py_format=False):
        vi = sys.version_info
        if py_format:
            return "py{0}{1}".format(vi.major, vi.minor)
        else:
            return "{0}.{1}".format(vi.major, vi.minor)

    def ersilia_version(self):
        ver = __ersilia_version__.split(".")
        ver = "{0}.{1}.{2}".format(ver[0], ver[1], ver[2].split("+")[0])
        return ver

    def ersilia_version_with_py(self):
        ver = self.ersilia_version()
        ver = "{0}-{1}".format(ver, self.python_version(py_format=True))
        return ver

    def ersilia_version_from_path(self, path):
        static_version_file = "_static_version.py"
        fn = os.path.join(path, "ersilia", static_version_file)
        if not os.path.exists(fn):
            fn = os.path.join(path, static_version_file)
        if not os.path.exists(fn):
            return None
        with open(fn, "r") as f:
            text = f.read()
            ver = text.split('"')[1]
        return ver

    def bentoml_version(self):
        return __bentoml_version__

    def server_docker_name(self, tag=None, as_tuple=False):
        if tag is None:
            tag = "{0}-{1}".format(
                self.ersilia_version(), self.python_version(py_format=True)
            )
        org = self.cfg.EXT.DOCKERHUB_ORG
        img = self.cfg.ENV.DOCKER.SERVER_BASE_IMAGE
        if as_tuple:
            return org, img, tag
        else:
            name = "{0}/{1}:{2}".format(img, org, tag)
            return name

    def base_conda_name(self, org, tag):
        if tag is None:
            tag = self.ersilia_version_with_py()
        env = self.cfg.ENV.CONDA.EOS_BASE_ENV
        name = "{0}-{1}-{2}".format(env, org, tag)
        return name

    @staticmethod
    def reformat_py(v):
        if len(v) < 4:
            raise Exception
        return "{0}.{1}".format(v[2], v[3:])
