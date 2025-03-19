import importlib.util
import os
import sys
from threading import Lock
from typing import Optional

from ...default import EOS
from ...tools.bentoml.exceptions import BentoMLException
from ...utils.logging import logger
from ...utils.terminal import run_command
from .setuptools_req import verify_setuptools


class BentoMLRequirement(object):
    """Handles installation and version checking of BentoML for Ersilia."""

    _lock = Lock()

    def __init__(self):
        verify_setuptools()
        self.logger = logger

    def is_installed(self) -> bool:
        """
        Checks if BentoML is installed.
        """
        module_name = "bentoml"
        if importlib.util.find_spec(module_name) is not None:
            return True
        else:
            self.logger.debug("BentoML is not installed")
            return False

    def _get_bentoml_version(self) -> Optional[str]:
        version_file = os.path.join(EOS, "bentomlversion.txt")

        if not os.path.exists(version_file):
            try:
                self.logger.debug("Importing BentoML module to resolve the version...")
                bentoml = importlib.import_module("bentoml")
                version = bentoml.__version__
                with open(version_file, "w") as f:
                    f.write("bentoml, version {0}".format(version))
                print(f"BentoML version: {version}")
            except:
                self.logger.warning("BentoML Python module could not be imported")
                version = None

            if version is None:
                cmd = "bentoml --version > {0}".format(version_file)
                run_command(cmd)

        with open(version_file, "r") as f:
            version_str = f.read().strip()

        if not version_str.startswith("bentoml, version"):
            raise Exception(
                "BentoML version was not stored correctly. Please check the bentomlversion.txt file in the eos folder"
            )

        version_str = version_str.split("version")[-1].strip()
        return version_str

    def is_bentoml_ersilia_version(self) -> bool:
        """
        Checks if the installed BentoML version is the Ersilia version
        ."""
        version_str = self._get_bentoml_version()
        if not version_str:
            return False
        if "0.11.0" in version_str:
            return True
        return False

    def _cleanup_corrupted_bentoml(self) -> None:
        with self._lock:
            self.logger.info("Cleaning up corrupted BentoML installation...")

            if self.is_installed():
                cmd = f"{sys.executable} -m pip uninstall bentoml -y"
                run_command(cmd)
            self.install(retries=1)

    def install(self, retries: int = 3) -> None:
        """
        Installs the Ersilia version of BentoML with error handling.
        """
        with self._lock:
            if retries <= 0:
                self.logger.critical(
                    "Final installation attempt failed. Manual intervention required."
                )
                raise BentoMLException("Installation failed after multiple attempts.")

            try:
                if self.is_installed() and not self.is_bentoml_ersilia_version():
                    self.logger.info("Uninstalling incompatible BentoML...")

                    cmd = f"{sys.executable} -m pip uninstall bentoml -y"
                    run_command(cmd)
                cmd = f"{sys.executable} -m pip install -U git+https://github.com/ersilia-os/bentoml-ersilia.git"
                run_command(cmd)
                if not self.is_bentoml_ersilia_version():
                    raise BentoMLException("Installed version verification failed")

                self.logger.info("BentoML was installed and verified successfully.")

            except BentoMLException as e:
                self.logger.error(f"Install attempt failed: {str(e)}")
                self.logger.warning(f"Attempts remaining: {retries - 1}")
                try:
                    self._cleanup_corrupted_bentoml()
                except Exception as cleanup_error:
                    self.logger.critical(f"Cleanup failed: {str(cleanup_error)}")
                    raise BentoMLException(
                        "Aborting due to failed cleanup"
                    ) from cleanup_error
                self.install(retries=retries - 1)
