import os
import subprocess
import sys

from ...default import EOS


class BentoMLRequirement(object):
    """
    A class to handle the installation and version checking of BentoML for Ersilia.

    Methods
    -------
    is_installed()
        Checks if BentoML is installed.
    is_bentoml_ersilia_version()
        Checks if the installed BentoML version is the Ersilia version.
    install()
        Installs the Ersilia version of BentoML.
    """

    def __init__(self):
        pass

    def is_installed(self) -> bool:
        """
        Checks if BentoML is installed.

        Returns
        -------
        bool
            True if BentoML is installed, False otherwise.
        """
        try:
            import bentoml  # noqa: F401

            return True
        except ImportError:
            return False

    def is_bentoml_ersilia_version(self) -> bool:
        """
        Checks if the installed BentoML version is the Ersilia version.

        Returns
        -------
        bool
            True if the installed BentoML version is the Ersilia version, False otherwise.
        """
        if not self.is_installed():
            return False

        version_file = os.path.join(EOS, "bentomlversion.txt")

        if not os.path.exists(version_file):
            cmd = "bentoml --version > {0}".format(version_file)
            subprocess.Popen(cmd, shell=True).wait()

        with open(version_file, "r") as f:
            text = f.read()
        if "0.11.0" in text:
            return True
        else:
            return False

    def install(self) -> None:
        """
        Installs the Ersilia version of BentoML.

        This method installs the BentoML package from the Ersilia GitHub repository.

        Returns
        -------
        None
        """
        print("Installing bentoml (the ersilia version)")
        cmd = "{0} -m pip install -U git+https://github.com/ersilia-os/bentoml-ersilia.git".format(
            sys.executable
        )
        subprocess.Popen(cmd, shell=True).wait()
