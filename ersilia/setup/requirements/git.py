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
    def is_installed(self, raise_exception=False, install_if_necessary=False):
        check = run_command_check_output("gh")
        if "GitHub" in check:
            return True
        else:
            if raise_exception:
                if install_if_necessary:
                    self.install()
                else:
                    raise GithubCliSetupError
            else:
                return False

    def install(self):
        run_command("conda install -c conda-forge gh")


class GitLfsRequirement(object):
    def __init__(self):
        self.name = "git-lfs"

    @throw_ersilia_exception
    def is_installed(self, install_if_necessary=True):
        check = run_command_check_output("git-lfs")
        if check.startswith("git-lfs"):
            return True
        else:
            if install_if_necessary:
                self.install()
            else:
                raise GitLfsSetupError

    def activate(self):
        run_command("git-lfs install")

    def install(self):
        run_command("conda install -c conda-forge git-lfs")
