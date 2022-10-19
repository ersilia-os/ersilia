import subprocess as s

from ... import throw_ersilia_exception
from ...utils.exceptions_utils.setup_exceptions import GitLfsSetupError


class GithubCliRequirement(object):
    def __ini__(self):
        self.name = "gh"
        if not self.is_installed():
            self.install()

    def is_installed(self):
        pass

    def install(self):
        pass


class GitLfsRequirement(object):
    def __init__(self):
        self.is_installed()
        self.activate()

    @throw_ersilia_exception
    def is_installed(self):
        try:
            check = s.run(["git-lfsxx"], capture_output=True)
        except:
            raise GitLfsSetupError

    def activate(self):
        s.run(["git-lfs install"])
