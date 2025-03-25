import os
import shutil
import sys

import click

from .. import ErsiliaBase, check_install_status
from ..default import CONFIG_JSON, EOS
from ..setup.baseconda import SetupBaseConda
from .conda import SimpleConda
from .config import Checker
from .logging import make_temp_dir
from .terminal import run_command
from .versioning import Versioner

INSTALL_LOG_FILE = ".install.log"


class BaseInstaller(ErsiliaBase):
    """
    Base class for managing installations in Ersilia.

    Parameters
    ----------
    check_install_log : bool
        Whether to check the installation log.
    config_json : dict, optional
        Configuration settings in JSON format. Default is None.
    credentials_json : dict, optional
        Credentials settings in JSON format. Default is None.
    """

    def __init__(self, check_install_log, config_json, credentials_json):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.check_install_log = check_install_log
        self.log_file = os.path.join(EOS, INSTALL_LOG_FILE)
        self.log = None
        self.read_log()
        self.versions = Versioner()
        checker = Checker()
        checker._config()
        checker._credentials()
        checker._package_path()
        self.development_path = checker.development_path

    def write_log(self):
        """
        Write the installation log to a file.
        """
        if self.log is None:
            return
        with open(self.log_file, "w") as f:
            for l in sorted(self.log):
                f.write(l + "\n")

    def update_log(self, task):
        """
        Update the installation log with a new task.

        Parameters
        ----------
        task : str
            The task to add to the log.
        """
        if self.log is None:
            self.log = {task}
        self.log.update([task])
        self.write_log()

    def remove_from_log(self, task):
        """
        Remove a task from the installation log.

        Parameters
        ----------
        task : str
            The task to remove from the log.
        """
        if self.log is not None:
            if task in self.log:
                self.log.remove(task)
                self.write_log()

    def read_log(self):
        """
        Read the installation log from a file.
        """
        if not os.path.exists(self.log_file):
            return
        with open(self.log_file, "r") as f:
            self.log = []
            for l in f:
                self.log += [l.rstrip()]
        self.log = set(self.log)

    def remove_log(self):
        """
        Remove the installation log file.
        """
        if os.path.exists(self.log_file):
            os.remove(self.log_file)

    def _is_done(self, name):
        if not self.check_install_log:
            return False
        if self.log is None:
            pass
        else:
            if name in self.log:
                return True
            else:
                pass
        self.update_log(name)
        return False

    @staticmethod
    def _is_tool(name):
        return shutil.which(name) is not None


class Installer(BaseInstaller):
    """
    A class to manage the installation of dependencies for Ersilia.

    Parameters
    ----------
    check_install_log : bool, optional
        Whether to check the installation log. Default is True.
    config_json : dict, optional
        Configuration settings in JSON format. Default is None.
    credentials_json : dict, optional
        Credentials settings in JSON format. Default is None.
    """

    def __init__(self, check_install_log=True, config_json=None, credentials_json=None):
        BaseInstaller.__init__(
            self,
            check_install_log=check_install_log,
            config_json=config_json,
            credentials_json=credentials_json,
        )

    def profile(self):
        """
        Set up the 'ersilia' CLI in the user profile.
        """
        if self._is_done("profile"):
            return
        from ..default import bashrc_cli_snippet

        click.echo(">> Setting up 'ersilia' CLI in user profile")
        bashrc_cli_snippet()

    def conda(self):
        """
        Check if Conda is installed.
        """
        if self._is_done("conda"):
            return
        if self._is_tool("conda"):
            return
        click.echo("Conda needs to be installed")
        sys.exit(1)

    def git(self):
        """
        Check if Git is installed.
        """
        if self._is_done("git"):
            return
        if self._is_tool("git"):
            return
        click.echo("Git needs to be installed")
        sys.exit(1)

    def config(self):
        """
        Set up the configuration file.
        """
        if self._is_done("config"):
            return
        if os.path.exists(os.path.join(EOS, CONFIG_JSON)):
            return
        click.echo(">> Setting up Config file")
        checker = Checker()
        checker.config()

    def _clone_repo(self, path):
        path_repo = os.path.join(path, self.cfg.HUB.PACKAGE)
        dev_path = self.development_path
        if dev_path is not None:
            shutil.copytree(dev_path, path_repo)
        else:
            from .download import GitHubDownloader

            gd = GitHubDownloader(overwrite=True)
            gd.clone(self.cfg.HUB.ORG, self.cfg.HUB.PACKAGE, path_repo)
        return path_repo

    def base_conda(self):
        """
        Create a base Conda environment.
        """
        if self._is_done("base_conda"):
            return
        eos_base_env = self.versions.base_conda_name()
        sc = SimpleConda()
        if sc.exists(eos_base_env):
            return
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_repo = self._clone_repo(tmp_folder)
        tmp_script = os.path.join(tmp_folder, "script.sh")
        tmp_python_script = os.path.join(tmp_folder, "base_installer.py")
        is_base = sc.is_base()
        if not is_base:
            bash_script = """
            source {0}/etc/profile.d/conda.sh
            conda deactivate
            """.format(sc.conda_prefix(False))
        else:
            bash_script = ""
        bash_script += """
        source {0}/etc/profile.d/conda.sh
        """.format(sc.conda_prefix(True))
        bc = SetupBaseConda()
        python_version = self.versions.python_version()
        python_version = bc.find_closest_python_version(python_version)
        bash_script += """
        cd {0}
        conda create -n {1} python={2} -y
        conda activate {1}
        pip install -e .
        python {3}
        conda deactivate
        """.format(tmp_repo, eos_base_env, python_version, tmp_python_script)
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        python_script = """
        from ersilia.utils.installers import base_installer
        base_installer(ignore_status=True)
        """
        with open(tmp_python_script, "w") as f:
            lines = python_script.split("\n")
            for l in lines:
                f.write(l[8:] + "\n")
        click.echo(">> Creating a Base Conda environment {0}".format(eos_base_env))
        run_command("bash {0}".format(tmp_script))

    def base_conda_slim(self):
        """
        Create a slim base Conda environment.
        """
        if self._is_done("base_conda_slim"):
            return
        # TODO

    def server_docker(self):
        """
        Build the Docker server image.
        """
        if self._is_done("server_docker"):
            return
        from .docker import SimpleDocker

        docker = SimpleDocker()
        org, img, tag = self.versions.server_docker_name(as_tuple=True)
        if docker.exists(org, img, tag):
            return
        # get a copy of the repository in a temporary directory
        tmp_dir = make_temp_dir(prefix="ersilia-")
        tmp_repo = self._clone_repo(tmp_dir)
        # write the dockerfile
        dockerfile = """
        FROM bentoml/model-server:{0}-{1}
        MAINTAINER ersilia

        ENV LC_ALL=C.UTF-8
        ENV LANG=C.UTF-8

        WORKDIR {2}

        COPY . .

        RUN conda --version

        RUN pip install .
        RUN ersilia setup --base
        """.format(
            self.versions.bentoml_version(),
            self.versions.python_version(py_format=True),
            self.cfg.ENV.DOCKER.IMAGE_WORKDIR,
        )
        path = os.path.join(tmp_repo, "Dockerfile")
        with open(path, "w") as f:
            lines = dockerfile.split("\n")
            lines = lines[1:-1]
            for l in lines:
                f.write(l[8:] + "\n")
        click.echo(
            ">> Building docker server image {0}".format(
                self.versions.server_docker_name(as_tuple=False)
            )
        )
        docker.build(path=tmp_repo, org=org, img=img, tag=tag)

    def server_docker_slim(self):
        """
        Build the slim Docker server image.
        """
        if self._is_done("server_docker_slim"):
            return
        # TODO


class Uninstaller(BaseInstaller):
    """
    A class to manage the uninstallation of dependencies for Ersilia.

    Parameters
    ----------
    check_install_log : bool, optional
        Whether to check the installation log. Default is True.
    config_json : dict, optional
        Configuration settings in JSON format. Default is None.
    credentials_json : dict, optional
        Credentials settings in JSON format. Default is None.
    """

    def __init__(self, check_install_log=True, config_json=None, credentials_json=None):
        BaseInstaller.__init__(
            self,
            check_install_log=check_install_log,
            config_json=config_json,
            credentials_json=credentials_json,
        )

    def base_conda(self):
        """
        Delete the base Conda environment.
        """
        self.remove_from_log("base_conda")
        eos_base_env = self.versions.base_conda_name()
        sc = SimpleConda()
        sc.delete(eos_base_env)

    def server_docker(self):
        """
        Delete the Docker server image.
        """
        self.remove_from_log("server_docker")
        from .docker import SimpleDocker

        docker = SimpleDocker()
        org, img, tag = self.versions.server_docker_name(as_tuple=True)
        if docker.exists(org, img, tag):
            docker.delete(org, img, tag)


def base_installer(ignore_status=False):
    """
    Perform a bare minimum installation of dependencies.

    This function is mainly used to create a base environment for the models.

    Parameters
    ----------
    ignore_status : bool, optional
        Whether to ignore the current installation status. Default is False.

    Returns
    -------
    None
    """
    status = check_install_status()
    if status["status"] is None or ignore_status:
        ins = Installer(check_install_log=False)
        ins.config()
        with open(status["install_status_file"], "w") as f:
            f.write("base")


def full_installer(ignore_status=False):
    """
    Perform a full installation of dependencies for Ersilia.

    Parameters
    ----------
    ignore_status : bool, optional
        Whether to ignore the current installation status. Default is False.
    """
    status = check_install_status()
    if status["status"] != "full" or ignore_status:
        ins = Installer()
        ins.profile()
        ins.conda()
        ins.git()
        ins.config()
        ins.base_conda()
        ins.server_docker()
        with open(status["install_status_file"], "w") as f:
            f.write("full")
