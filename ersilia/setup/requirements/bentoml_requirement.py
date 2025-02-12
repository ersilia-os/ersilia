import sys
from threading import Lock
from typing import Optional

from packaging import version
from packaging.version import InvalidVersion

from ...tools.bentoml.exceptions import BentoMLException
from ...utils.logging import logger
from ...utils.terminal import run_command


class BentoMLRequirement(object):
    """Handles installation and version checking of BentoML for Ersilia."""

    _lock = Lock()

    def __init__(self):
        self.logger = logger

    def is_installed(self) -> bool:
        """Checks if BentoML is installed."""
        try:
            import bentoml  # noqa: F401

            return True
        except ImportError:
            self.logger.debug("BentoML is not installed")
            return False

    def _get_bentoml_version(self) -> Optional[str]:
        """Get BentoML version using run_command"""
        result = run_command([sys.executable, "-m", "bentoml", "--version"])

        if result.returncode != 0:
            self.logger.error(f"BentoML version check failed: {result.stderr}")
            return None

        # Extract version from "python -m bentoml, version 0.11.0"
        version_str = result.stdout.split("version")[-1].strip()
        return version_str

        # try:
        #    return version.parse(version_str).public
        # except InvalidVersion:
        #    self.logger.error(f"Invalid BentoML version detected: {version_str}")
        #    return None

    def is_bentoml_ersilia_version(self) -> bool:
        """Checks if the installed BentoML version is the Ersilia version."""
        version_str = self._get_bentoml_version()
        if not version_str:
            return False
        try:
            return version.parse(version_str) == version.parse("0.11.0")
        except InvalidVersion:
            self.logger.error(f"Invalid BentoML version detected: {version_str}")
            return False

    def _cleanup_corrupted_bentoml(self) -> None:
        """Forcefully uninstall BentoML and reinstall the correct version."""
        with self._lock:
            self.logger.info("Cleaning up corrupted BentoML installation...")

            if self.is_installed():
                result = run_command(
                    [sys.executable, "-m", "pip", "uninstall", "bentoml", "-y"]
                )
                if result.returncode != 0:
                    raise BentoMLException(f"Force uninstall failed: {result.stderr}")

            # ✅ Ensure reinstallation after cleanup
            self.logger.info("Reinstalling Ersilia-compatible BentoML after cleanup...")
            self.install(retries=1)  # Only allow 1 retry to avoid infinite loops

    def install(self, retries: int = 3) -> None:
        """Installs the Ersilia version of BentoML with error handling."""
        with self._lock:
            if retries <= 0:
                self.logger.critical(
                    "Final installation attempt failed. Manual intervention required."
                )
                raise BentoMLException("Installation failed after multiple attempts.")

            try:
                # 1. Uninstall if wrong version exists
                if self.is_installed() and not self.is_bentoml_ersilia_version():
                    self.logger.info("Uninstalling incompatible BentoML...")
                    uninstall_result = run_command(
                        [sys.executable, "-m", "pip", "uninstall", "bentoml", "-y"]
                    )
                    if uninstall_result.returncode != 0:
                        raise BentoMLException(
                            f"Uninstall failed: {uninstall_result.stderr}"
                        )

                # 2. Install specific version
                self.logger.info("Installing Ersilia-compatible BentoML...")
                install_result = run_command(
                    [
                        sys.executable,
                        "-m",
                        "pip",
                        "install",
                        "-U",
                        "git+https://github.com/ersilia-os/bentoml-ersilia.git",
                    ]
                )
                if install_result.returncode != 0:
                    raise BentoMLException(f"Install failed: {install_result.stderr}")

                # 3. Post-install verification
                if not self.is_bentoml_ersilia_version():
                    raise BentoMLException("Installed version verification failed")

                self.logger.info("BentoML was installed and verified successfully.")

            except BentoMLException as e:
                self.logger.error(f"Install attempt failed: {str(e)}")
                self.logger.warning(f"Attempts remaining: {retries-1}")

                try:
                    self._cleanup_corrupted_bentoml()
                except Exception as cleanup_error:
                    self.logger.critical(f"Cleanup failed: {str(cleanup_error)}")
                    raise BentoMLException(
                        "Aborting due to failed cleanup"
                    ) from cleanup_error

                # Recursive retry with counter
                self.install(retries=retries - 1)
