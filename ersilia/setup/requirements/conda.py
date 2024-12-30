from ...utils.terminal import run_command_check_output


class CondaRequirement(object):
    """
    A class to handle the installation and checking of the Conda package manager.

    Methods
    -------
    is_installed()
        Checks if Conda is installed.
    install()
        Placeholder for installing Conda.
    """

    def __init__(self):
        self.name = "conda"

    def is_installed(self) -> bool:
        """
        Checks if Conda is installed.

        Returns
        -------
        bool
            True if Conda is installed, False otherwise.
        """
        cmd = "command -v {0}".format(self.name)
        output = run_command_check_output(cmd)
        if self.name in output:
            return True
        else:
            return False

    def install(self) -> None:
        """
        Placeholder for installing Conda.

        Returns
        -------
        None
        """
        pass
        # TODO
