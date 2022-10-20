from ... import throw_ersilia_exception
from ...utils.exceptions_utils.setup_exceptions import (
    GitLfsSetupError,
    GithubCliSetupError,
)
from ...utils.terminal import run_command, run_command_check_output


class GithubCliRequirement(object):
    def __ini__(self):
        self.name = "gh"

    @throw_ersilia_exception
    def is_installed(self, raise_exception):
        check = run_command_check_output("gh")
        if "GitHub" in check:
            return True
        else:
            if raise_exception:
                raise GithubCliSetupError
            else:
                return False


class GitLfsRequirement(object):
    def __init__(self):
        self.name = "git-lfs"

    @throw_ersilia_exception
    def is_installed(self):
        check = run_command_check_output("git-lfs")
        if check.startswith("git-lfs"):
            return True
        else:
            raise GitLfsSetupError

    def activate(self):
        run_command("git-lfs install")
