import os
import shutil

from .. import ErsiliaBase, logger, throw_ersilia_exception
from ..utils.exceptions_utils.fetch_exceptions import (
    ModelPackageInstallError,
    VirtualEnvironmentSetupError,
)
from ..utils.logging import make_temp_dir
from .terminal import run_command


class SimpleVenv(ErsiliaBase):
    """
    A class to manage virtual environments for Ersilia.

    Parameters
    ----------
    root : str
        The root directory for the virtual environments.

    Methods
    -------
    create(environment)
        Create a new virtual environment.
    delete(environment)
        Delete a virtual environment.
    run_commandlines(environment, commandlines)
        Run command lines in a virtual environment.
    """

    def __init__(self, root):
        ErsiliaBase.__init__(self, config_json=None, credentials_json=None)
        self.root = os.path.abspath(root)
        self.logger.debug("Setting virtual environment at {0}".format(self.root))

    def _get_path(self, environment):
        return os.path.join(self.root, environment)

    def exists(self, environment):
        """
        Check if a virtual environment exists.

        Parameters
        ----------
        environment : str
            The name of the virtual environment.

        Returns
        -------
        bool
            True if the virtual environment exists, False otherwise.
        """
        if os.path.exists(self._get_path(environment)):
            return True
        else:
            return False

    @throw_ersilia_exception()
    def create(self, environment):
        """
        Create a new virtual environment.

        Parameters
        ----------
        environment : str
            The name of the virtual environment.

        Returns
        -------
        None

        Raises
        ------
        VirtualEnvironmentSetupError
            If the virtual environment setup fails.
        """
        path = self._get_path(environment)
        if self.exists(path):
            return
        run_command("python -m venv {0} --symlinks --clear".format(path))
        if not self.exists(environment):
            raise VirtualEnvironmentSetupError(environment)

    def delete(self, environment):
        """
        Delete a virtual environment.

        Parameters
        ----------
        environment : str
            The name of the virtual environment.

        Returns
        -------
        None
        """
        path = self._get_path(environment)
        if not self.exists(path):
            return
        shutil.rmtree(path)

    @throw_ersilia_exception()
    def run_commandlines(self, environment, commandlines):
        """
        Run command lines in a virtual environment.

        Parameters
        ----------
        environment : str
            The name of the virtual environment.
        commandlines : str
            The command lines to run.

        Returns
        -------
        None

        Raises
        ------
        ModelPackageInstallError
            If the command lines execution fails.
        """
        if not self.exists(environment):
            raise Exception("{0} environment does not exist".format(environment))
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_script = os.path.join(tmp_folder, "script.sh")
        tmp_log = os.path.join(tmp_folder, "installs.log")  # new
        with open(tmp_script, "w") as f:
            f.write("cd {0}{1}".format(self.root, os.linesep))
            f.write("source {0}/bin/activate{1}".format(environment, os.linesep))
            f.write("{0}{1}".format(commandlines, os.linesep))
            f.write("deactivate")

        logger.debug("Running {0}".format(tmp_script))  #

        # run_command("bash {0}".format(tmp_script))
        run_command("bash {0} 2>&1 | tee -a {1}".format(tmp_script, tmp_log))  #
        with open(tmp_log, "r") as f:
            log_file = f.read()
        logger.debug(log_file)
        if "ERROR" in log_file:
            logger.debug("venv ModelPackageInstallError")
            raise ModelPackageInstallError(tmp_script)
        logger.debug("Activation done")
