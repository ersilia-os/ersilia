import os
import yaml

from . import BasePack
from .....utils.terminal import run_command
from .....db.environments.localdb import EnvironmentDb
from .....utils.conda import SimpleConda
from .....utils.docker import SimpleDockerfileParser

from ..... import throw_ersilia_exception
from .....utils.exceptions_utils.fetch_exceptions import CondaEnvironmentExistsError

from ..... import EOS


class SystemPack(BasePack):
    def __init__(self, model_id, config_json):
        BasePack.__init__(self, model_id, config_json)
        self.logger.debug("Initializing system packer")

    def _run(self):
        self.logger.debug("Packing model with system installation")
        dest_dir = self._model_path(self.model_id)
        bundle_dir = os.path.join(EOS, "repository")
        cmd = "ersilia_model_pack --repo_path {0} --bundles_repo_path {1}".format(
            dest_dir, bundle_dir
        )
        self.logger.debug("Running packing command {0}".format(cmd))
        run_command(cmd)
        self._symlinks()

    def run(self):
        self._run()


class CondaPack(BasePack):
    def __init__(self, model_id, config_json):
        BasePack.__init__(self, model_id, config_json)
        self.conda = SimpleConda()
        self.logger.debug("Initializing conda packer")

    def _decide_python_version(self):
        install_yml_file = os.path.join(self._model_path(self.model_id), "install.yml")
        if os.path.exists(install_yml_file):
            self.logger.debug("Reading python version from install.yml")
            with open(install_yml_file, "r") as f:
                data = yaml.safe_load(f)
                return data["python"]
        dockerfile_file = os.path.join(self._model_path(self.model_id), "Dockerfile")
        if os.path.exists(dockerfile_file):
            self.logger.debug("Reading python version from Dockerfile")
            pyver = (
                SimpleDockerfileParser(dockerfile_file).get_baseimage().split("-py")[-1]
            )
            pyver = pyver[0] + "." + pyver[1:]
            return pyver
        self.logger.error("Could not find python version in install.yml or Dockerfile")
        return None

    def _setup(self):
        self.logger.debug("Setting up")
        model_id = self.model_id
        env = model_id
        python_version = self._decide_python_version()
        self.logger.debug("Conda environment {0}".format(env))
        if not self.conda.exists(env):
            self.logger.debug("Environment {0} does not exist".format(env))
            self.conda.create(environment=env, python_version=python_version)
            self.logger.debug("Creating base conda environment")
            commandlines = [
                "python -m pip install git+https://github.com/ersilia-os/ersilia-pack.git"
            ]
            self.conda.run_commandlines(environment=env, commandlines=commandlines)
            self.logger.debug(
                "Storing Conda environment in the local environment database"
            )
            db = EnvironmentDb(config_json=self.config_json)
            db.table = "conda"
            db.insert(model_id=model_id, env=env)
            self.logger.debug("Done with the Conda setup")
        else:
            self.logger.debug("Environment {0} does exist".format(env))

        if env is not None:
            if not self.conda.exists(env):
                raise CondaEnvironmentExistsError(env)
        return env

    def _run(self):
        env = self._setup()
        self.logger.debug("Using environment {0}".format(env))
        dest_dir = self._model_path(self.model_id)
        bundle_dir = os.path.join(EOS, "repository")
        cmd = "ersilia_model_pack --repo_path {0} --bundles_repo_path {1} --conda_env_name {2}".format(
            dest_dir, bundle_dir, env
        )
        self.logger.debug("Running command: {0}".format(cmd))
        self.conda.run_commandlines(environment=env, commandlines=cmd)
        self.logger.debug(
            "Packing command successfully run inside {0} conda environment".format(env)
        )
        self._symlinks()
        self.logger.debug("Symlinks created")

    @throw_ersilia_exception
    def run(self):
        self.logger.debug("Packing model with Conda")
        self._run()


def get_runner(pack_mode):
    if pack_mode == "system":
        return SystemPack
    if pack_mode == "conda":
        return CondaPack
