import importlib

from ...utils.terminal import run_command


class IsauraRequirement(object):
    """
    A class to handle the installation of the Isaura library.

    Methods
    -------
    install()
        Installs the Isaura library.
    """

    def __init__(self):
        self.name = "isaura"
        try:
            importlib.import_module(self.name)
        except:
            self.install()

    def install(self) -> None:
        """
        Installs the Isaura library.

        This method installs the Isaura package using a placeholder command.

        Returns
        -------
        None
        """
        run_command()
