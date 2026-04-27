import importlib
import sys

from ... import throw_ersilia_exception
from ...utils.exceptions_utils.setup_exceptions import (
    GithubCliSetupError,
    GitLfsSetupError,
)
from ...utils.terminal import run_command, run_command_check_output


class GithubCliRequirement(object):
    """
    A class to handle the checking and installation of the GitHub CLI.

    Methods
    -------
    is_installed(raise_exception=False, install_if_necessary=False)
        Checks if the GitHub CLI is installed.
    install()
        Installs the GitHub CLI.
    """

    def __init__(self):
        self.name = "gh"

    @throw_ersilia_exception()
    def is_installed(self, raise_exception=False, install_if_necessary=False) -> bool:
        """
        Checks if the GitHub CLI is installed.

        Parameters
        ----------
        raise_exception : bool, optional
            Whether to raise an exception if the GitHub CLI is not installed (default is False).
        install_if_necessary : bool, optional
            Whether to install the GitHub CLI if it is not installed (default is False).

        Returns
        -------
        bool
            True if the GitHub CLI is installed, False otherwise.
        """
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

    def install(self) -> None:
        """
        Installs the GitHub CLI.

        This method installs the GitHub CLI using conda.

        Returns
        -------
        None
        """
        run_command("conda install -c conda-forge gh")


class GitLfsRequirement(object):
    """
    A class to handle the checking and installation of Git LFS.

    Methods
    -------
    is_installed(install_if_necessary=True)
        Checks if Git LFS is installed.
    activate()
        Activates Git LFS.
    install()
        Installs Git LFS.
    """

    def __init__(self):
        self.name = "git-lfs"

    @throw_ersilia_exception()
    def is_installed(self, install_if_necessary=True) -> bool:
        """
        Checks if Git LFS is installed.

        Parameters
        ----------
        install_if_necessary : bool, optional
            Whether to install Git LFS if it is not installed (default is True).

        Returns
        -------
        bool
            True if Git LFS is installed, False otherwise.
        """
        check = run_command_check_output("git-lfs")
        if check.startswith("git-lfs"):
            return True
        else:
            if install_if_necessary:
                self.install()
            else:
                raise GitLfsSetupError

    def activate(self) -> None:
        """
        Activates Git LFS.

        This method runs the Git LFS install command.

        Returns
        -------
        None
        """
        run_command("git-lfs install")

    def install(self) -> None:
        """
        Installs Git LFS.

        This method installs Git LFS using conda.

        Returns
        -------
        None
        """
        run_command("conda install -c conda-forge git-lfs")


class EosvcRequirement(object):
    """
    A class to handle checking and installing eosvc in the active Python environment.
    """

    def __init__(self):
        self.name = "eosvc"
        self.install_target = "git+https://github.com/ersilia-os/eosvc.git"

    def is_installed(self) -> bool:
        """
        Check whether eosvc is importable in the active environment.
        """
        try:
            importlib.import_module(self.name)
            return True
        except ModuleNotFoundError:
            return False

    def install(self) -> bool:
        """
        Install eosvc into the active Python environment.
        """
        cmd = [
            sys.executable,
            "-m",
            "pip",
            "install",
            self.install_target,
        ]
        result = run_command(cmd)
        return result.returncode == 0

    def ensure_installed(self) -> bool:
        """
        Ensure eosvc is available in the active Python environment.
        """
        if self.is_installed():
            return True
        if not self.install():
            return False
        return self.is_installed()
