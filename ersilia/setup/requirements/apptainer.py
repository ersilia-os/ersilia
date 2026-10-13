import platform
import shutil

from ...utils.exceptions_utils.cli_exceptions import (
    ApptainerNotInstalledError,
    ApptainerNotLinuxError,
)


class ApptainerRequirement(object):
    """
    Checks that Apptainer images can be run on this machine.

    Apptainer (formerly Singularity) runs on Linux only. Both the ``apptainer``
    and the older ``singularity`` command are accepted; they take the same
    arguments.

    Methods
    -------
    is_linux()
        Checks if this machine runs Linux.
    binary()
        The Apptainer command available, if any.
    is_installed()
        Checks if Apptainer is installed.
    check()
        Raises a clear error if Apptainer images cannot be run here.
    """

    BINARIES = ("apptainer", "singularity")

    def __init__(self):
        self.name = "apptainer"

    def is_linux(self) -> bool:
        """
        Checks if this machine runs Linux.

        Returns
        -------
        bool
            True on Linux, False otherwise.
        """
        return platform.system() == "Linux"

    def binary(self):
        """
        The Apptainer command available on the PATH.

        Returns
        -------
        str or None
            ``"apptainer"`` or ``"singularity"``, or None if neither is installed.
        """
        for name in self.BINARIES:
            if shutil.which(name):
                return name
        return None

    def is_installed(self) -> bool:
        """
        Checks if Apptainer is installed.

        Returns
        -------
        bool
            True if the ``apptainer`` or ``singularity`` command is available.
        """
        return self.binary() is not None

    def check(self):
        """
        Raises a clear error if Apptainer images cannot be run on this machine.

        Returns
        -------
        str
            The Apptainer command to use.

        Raises
        ------
        ApptainerNotLinuxError
            If this machine does not run Linux.
        ApptainerNotInstalledError
            If Apptainer is not installed.
        """
        if not self.is_linux():
            system = platform.system()
            raise ApptainerNotLinuxError(
                {"Darwin": "macOS"}.get(system, system or "this system")
            )
        binary = self.binary()
        if binary is None:
            raise ApptainerNotInstalledError()
        return binary
