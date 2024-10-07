import os
import shutil

from .conda import SimpleConda
from .docker import SimpleDocker
from ..default import EOS, BENTOML_PATH
from .logging import logger


class Uninstaller(object):
    def __init__(self):
        self.docker_cleaner = SimpleDocker()

    def _directories(self):
        dirs_to_remove = [EOS, BENTOML_PATH]
        for dir in dirs_to_remove:
            if os.path.exists(dir):
                try:
                    shutil.rmtree(dir)
                except Exception as e:
                    logger.error(f"Failed to remove directory {dir}")

    def _conda(self):
        sc = SimpleConda()
        for env in sc._env_list():
            if env.startswith("#"):
                continue
            if not env.startswith("eos"):
                continue
            env = env.split(" ")[0]
            if len(env.split("-")[0]) == 7:
                try:
                    logger.info(f"Removing conda environment: {env}")
                    sc.delete(env)
                except Exception as e:
                    logger.error(f"Failed to remove conda environment {env}: {e}")

    def uninstall(self):
        """Main uninstallation method"""

        try:
            logger.info("Starting Ersillia uninstallation...")

            self.docker_cleaner.cleanup_ersilia_images()
            self._conda()
            self._directories()
            logger.info("Ersilia uninstallation completed")

        except Exception as e:
            logger.error(f"Uninstallation failed: {e}")
            raise
