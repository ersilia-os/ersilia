from ....setup.requirements.eospath import EosHomePathRequirement
from ....setup.requirements.git import GitLfsRequirement, GithubCliRequirement
from ....setup.requirements.ping import PingRequirement
from ....setup.requirements.conda import CondaRequirement
from . import BaseAction


class SetupChecker(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )

    def _gh_cli(self):
        req = GithubCliRequirement()
        if req.is_installed(raise_exception=False):
            self.logger.debug("GitHub CLI is installed")
        else:
            self.logger.info(
                "GitHub CLI is not installed. Ersilia can work without it, but we highy recommend that you install this tool."
            )

    def _git_lfs(self):
        req = GitLfsRequirement()
        req.is_installed()
        self.logger.debug("Git LFS is installed")
        req.activate()
        self.logger.debug("Git LFS has been activated")

    def _ping(self):
        req = PingRequirement()
        req.is_connected()
        self.logger.debug("Connected to the internet")

    def _conda(self):
        req = CondaRequirement()
        req.is_installed()
        self.logger.debug("Conda is installed")

    def _eos_home_path(self):
        req = EosHomePathRequirement()
        req.eos_home_path_exists()
        self.logger.debug("EOS Home path exists")

    def check(self):
        self._gh_cli()
        self._git_lfs()
        self._ping()
        self._conda()
        self._eos_home_path()
