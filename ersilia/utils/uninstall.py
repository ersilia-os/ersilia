import os
import shutil
import subprocess

from ..default import BENTOML_PATH, EOS
from .conda import SimpleConda
from .docker import SimpleDocker
from .logging import logger


class Uninstaller(object):
    """
    A class to manage the uninstallation of Ersilia and its dependencies.

    Methods
    -------
    uninstall()
        Main uninstallation method.
    """

    def __init__(self):
        self.docker_cleaner = SimpleDocker()

    def _uninstall_ersilia_package(self):
        try:
            logger.info("Uninstalling Ersilia package...")
            subprocess.run(["pip", "uninstall", "-y", "ersilia"], check=True)
            logger.info("Ersilia package uninstalled successfully.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to uninstall Ersilia package: {e}")

    def _directories(self):
        dirs_to_remove = [EOS, BENTOML_PATH]
        for dir in dirs_to_remove:
            if os.path.exists(dir):
                try:
                    logger.info(f"Removing directory: {dir}...")
                    shutil.rmtree(dir)
                    logger.info(f"Directory {dir} removed successfully.")
                except Exception as e:
                    logger.error(f"Failed to remove directory {dir}: {e}")

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

        env_name = "ersilia"

        try:
            logger.info(f"Removing Conda environment: {env_name}...")
            sc.delete(env_name)
            logger.info(f"Conda environment {env_name} removed successfully.")
        except Exception as e:
            logger.error(f"Failed to remove Conda environment {env_name}: {e}")

    def uninstall(self):
        """
        Main uninstallation method.
        """
        try:
            logger.info("Starting Ersilia uninstallation...")

            self.docker_cleaner.cleanup_ersilia_images()
            self._uninstall_ersilia_package()
            self._conda()
            self._directories()

            logger.info("Ersilia uninstallation completed")
        except Exception as e:
            logger.error(f"Uninstallation failed: {e}")
            raise
