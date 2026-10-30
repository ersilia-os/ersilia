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
        self.failed = []

    def _step(self, label, func):
        # The spinner prints the step as done or failed; the error goes below it.
        try:
            spinner(label, func)
        except Exception as e:
            self.failed.append(label[0].lower() + label[1:])
            echo(str(e))

    def _uninstall_ersilia_package(self):
        def _run():
            subprocess.run(["pip", "uninstall", "-y", "ersilia"], check=True)

        self._step("Uninstalling the ersilia Python package", _run)

    def _directories(self):
        def _run():
            for dir in [EOS]:
                if os.path.exists(dir):
                    shutil.rmtree(dir)

        self._step(f"Removing the Ersilia folder ({EOS})", _run)

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

        self._step("Removing the model conda environments", _run)

    def uninstall(self):
        """
        Main uninstallation method.
        """
        self._step(
            "Removing the Ersilia Docker images",
            self.docker_cleaner.cleanup_ersilia_images,
        )
        self._uninstall_ersilia_package()
        self._conda()
        self._directories()

        if self.failed:
            echo(
                "Ersilia was not fully uninstalled. Failed steps: {0}.".format(
                    "; ".join(self.failed)
                ),
                fg="red",
            )
        else:
            echo("Ersilia uninstalled.", fg="green")
