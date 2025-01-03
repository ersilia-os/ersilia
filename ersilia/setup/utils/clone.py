import os
import shutil

from ... import ErsiliaBase
from ...utils.config import Checker
from ...utils.download import GitHubDownloader
from ...utils.versioning import Versioner


class ErsiliaCloner(ErsiliaBase):
    """
    A class to handle cloning of the Ersilia repository.

    Methods
    -------
    clone(path, version)
        Clones the Ersilia repository to the specified path and version.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        checker = Checker()
        checker._package_path()
        self.development_path = checker.get_development_path()

    def clone(self, path: str, version: str) -> str:
        """
        Clones the Ersilia repository to the specified path and version.

        Parameters
        ----------
        path : str
            The path where the repository should be cloned.
        version : str
            The version of the repository to clone.

        Returns
        -------
        str
            The path to the cloned repository.
        """
        path_repo = os.path.join(path, self.cfg.HUB.PACKAGE)
        if self.development_path is not None:
            path_version = Versioner().ersilia_version_from_path(self.development_path)
            if path_version is None:
                shutil.copytree(self.development_path, path_repo)
                return path_repo
            if path_version == version:
                shutil.copytree(self.development_path, path_repo)
                return path_repo
        gd = GitHubDownloader(overwrite=True)
        gd.clone(self.cfg.HUB.ORG, self.cfg.HUB.PACKAGE, path_repo)
        return path_repo
