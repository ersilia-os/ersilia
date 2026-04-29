import os
import shutil
import subprocess

from ..default import EOS
from .conda import SimpleConda
from .docker import SimpleDocker
from .echo import echo, spinner
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
        def _run():
            subprocess.run(["pip", "uninstall", "-y", "ersilia"], check=True)

        try:
            spinner("Uninstalling Ersilia pip package", _run)
        except subprocess.CalledProcessError as e:
            echo(f"Failed to uninstall Ersilia pip package: {e}", fg="red")

    def _directories(self):
        def _run():
            for dir in [EOS]:
                if os.path.exists(dir):
                    shutil.rmtree(dir)

        try:
            spinner(f"Removing EOS directory {EOS}", _run)
        except Exception as e:
            echo(f"Failed to remove EOS directory: {e}", fg="red")

    def _conda(self):
        sc = SimpleConda()

        def _run():
            for env in sc._env_list():
                if env.startswith("#"):
                    continue
                if not env.startswith("eos"):
                    continue
                env = env.split(" ")[0]
                if len(env.split("-")[0]) == 7:
                    try:
                        sc.delete(env)
                    except Exception as e:
                        logger.error(f"Failed to remove conda environment {env}: {e}")
            try:
                sc.delete("ersilia")
            except Exception as e:
                logger.error(f"Failed to remove conda environment ersilia: {e}")

        try:
            spinner("Removing model conda environments", _run)
        except Exception as e:
            echo(f"Failed to remove conda environments: {e}", fg="red")

    def uninstall(self):
        """
        Main uninstallation method.
        """
        try:
            spinner(
                "Removing Ersilia Docker images",
                self.docker_cleaner.cleanup_ersilia_images,
            )
        except Exception as e:
            echo(f"Failed to remove Docker images: {e}", fg="red")

        self._uninstall_ersilia_package()
        self._conda()
        self._directories()

        echo("Ersilia uninstalled successfully.", fg="green")
