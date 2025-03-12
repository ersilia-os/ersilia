from ....setup.requirements.conda import CondaRequirement
from ....setup.requirements.eospath import EosHomePathRequirement
from ....setup.requirements.git import GithubCliRequirement, GitLfsRequirement
from ....setup.requirements.ping import PingRequirement
from . import BaseAction


class SetupChecker(BaseAction):
    """
    Checks the setup requirements for the model to be installed and run.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    config_json : dict
        Configuration settings for the setup checker.

    Methods
    -------
    check()
        Checks all setup requirements.
    """

    def __init__(self, model_id: str, config_json: dict):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )

    def _gh_cli(self):
        req = GithubCliRequirement()
        if req.is_installed(raise_exception=False):
            self.logger.debug("GitHub CLI is installed")
        else:
            self.logger.info(
                "GitHub CLI is not installed. Ersilia can work without it, but we highly recommend that you install this tool."
            )

    def _git_lfs(self):
        self.logger.debug("Checking Git LFS")
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
        """
        Checks all setup requirements.
        """
        self._gh_cli()
        self._git_lfs()
        self._ping()
        self._conda()
        self._eos_home_path()
