import shutil
import os

from ... import ErsiliaBase
from ...utils.download import GitHubDownloader
from ...utils.config import Checker
from ...utils.versioning import Versioner


class ErsiliaCloner(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        checker = Checker()
        checker._package_path()
        self.development_path = checker.get_development_path()

    def clone(self, path, version):
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
