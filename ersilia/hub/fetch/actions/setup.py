from ersilia.setup.requirements.git import GitLfsRequirement, GithubCliRequirement
from . import BaseAction


class SetupChecker(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )

    def _gh_cli(self):
        req = GithubCliRequirement()
        req.is_installed()
        self.logger.debug("GitHub CLI is installed")

    def _git_lfs(self):
        req = GitLfsRequirement()
        req.is_installed()
        self.logger.debug("Git LFS is installed")
        req.activate()
        self.logger.debug("Git LFS has been activated")

    def check(self):
        self._gh_cli()
        self._git_lfs()
        #self._is_git
