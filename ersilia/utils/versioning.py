import os
import sys

try:
    from bentoml import __version__ as __bentoml_version__
except:
    __bentoml_version__ = None
from .. import ErsiliaBase
from .. import __version__ as __ersilia_version__


class Versioner(ErsiliaBase):
    """
    A class to manage versioning information for Ersilia.

    Methods
    -------
    python_version(py_format=False)
        Get the current Python version.
    ersilia_version()
        Get the current Ersilia version.
    ersilia_version_with_py()
        Get the current Ersilia version with Python version.
    ersilia_version_from_path(path)
        Get the Ersilia version from a given path.
    bentoml_version()
        Get the current BentoML version.
    server_docker_name(tag=None, as_tuple=False)
        Get the server Docker image name.
    base_conda_name(org, tag)
        Get the base Conda environment name.
    reformat_py(v)
        Reformat a Python version string.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def python_version(self, py_format=False):
        """
        Get the current Python version.

        Parameters
        ----------
        py_format : bool, optional
            Whether to return the version in Python format (e.g., 'py38'). Default is False.

        Returns
        -------
        str
            The current Python version.
        """
        vi = sys.version_info
        if py_format:
            return "py{0}{1}".format(vi.major, vi.minor)
        else:
            return "{0}.{1}".format(vi.major, vi.minor)

    def ersilia_version(self):
        """
        Get the current Ersilia version.

        Returns
        -------
        str
            The current Ersilia version.
        """
        ver = __ersilia_version__.split(".")
        ver = "{0}.{1}.{2}".format(ver[0], ver[1], ver[2].split("+")[0])
        return ver

    def ersilia_version_with_py(self):
        """
        Get the current Ersilia version with Python version.

        Returns
        -------
        str
            The current Ersilia version with Python version.
        """
        ver = self.ersilia_version()
        ver = "{0}-{1}".format(ver, self.python_version(py_format=True))
        return ver

    def ersilia_version_from_path(self, path):
        """
        Get the Ersilia version from a given path.

        Parameters
        ----------
        path : str
            The path to the Ersilia installation.

        Returns
        -------
        str or None
            The Ersilia version if found, otherwise None.
        """
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
        """
        Get the current BentoML version.

        Returns
        -------
        str or None
            The current BentoML version if available, otherwise None.
        """
        return __bentoml_version__

    def server_docker_name(self, tag=None, as_tuple=False):
        """
        Get the server Docker image name.

        Parameters
        ----------
        tag : str, optional
            The tag for the Docker image. Default is None.
        as_tuple : bool, optional
            Whether to return the name as a tuple. Default is False.

        Returns
        -------
        str or tuple
            The server Docker image name, or a tuple of (org, img, tag) if as_tuple is True.
        """
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
        """
        Get the base Conda environment name.

        Parameters
        ----------
        org : str
            The organization name.
        tag : str
            The tag for the Conda environment.

        Returns
        -------
        str
            The base Conda environment name.
        """
        if tag is None:
            tag = self.ersilia_version_with_py()
        env = self.cfg.ENV.CONDA.EOS_BASE_ENV
        name = "{0}-{1}-{2}".format(env, org, tag)
        return name

    @staticmethod
    def reformat_py(v):
        """
        Reformat a Python version string.

        Parameters
        ----------
        v : str
            The Python version string.

        Returns
        -------
        str
            The reformatted Python version string.

        Raises
        ------
        Exception
            If the version string is too short.
        """
        if len(v) < 4:
            raise Exception
        return "{0}.{1}".format(v[2], v[3:])
