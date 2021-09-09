import os
import tempfile

from . import BasePack
from ....utils.terminal import run_command
from ....db.environments.localdb import EnvironmentDb
from ....utils.venv import SimpleVenv
from ....utils.conda import SimpleConda
from ....utils.docker import SimpleDocker
from ....setup.baseconda import SetupBaseConda

from ....default import DEFAULT_VENV
from .. import PYTHON_INSTALLS


USE_CHECKSUM = False


class SystemPack(BasePack):
    def __init__(self, model_id, config_json):
        BasePack.__init__(self, model_id, config_json)
        self.logger.debug("Initializing system packer")

    def run(self):
        self.logger.debug("Packing model with system installation")
        folder = self._model_path(self.model_id)
        script_path = os.path.join(folder, self.cfg.HUB.PACK_SCRIPT)
        run_command("python {0}".format(script_path))
        self._bentoml_bundle_symlink()


class VenvPack(BasePack):
    def __init__(self, model_id, config_json):
        BasePack.__init__(self, model_id, config_json)
        self.logger.debug("Initializing virtualenv packer")

    def _setup(self):
        model_path = self._model_path(self.model_id)
        installs_file = os.path.join(model_path, PYTHON_INSTALLS)
        self.logger.debug("Reading python installs from {0}".format(installs_file))
        venv = SimpleVenv(model_path)
        venv.create(DEFAULT_VENV)
        with open(installs_file, "r") as f:
            commandlines = f.read()
        venv.run_commandlines(environment=DEFAULT_VENV, commandlines=commandlines)
        return venv

    def _run(self):
        venv = self._setup()
        pack_snippet = """
        python {0}
        """.format(
            self.cfg.HUB.PACK_SCRIPT
        )
        venv.run_commandlines(environment=DEFAULT_VENV, commandlines=pack_snippet)
        self._bentoml_bundle_symlink()

    def run(self):
        self.logger.debug("Packing model with VirtualEnv")
        self._write_python_installs()
        self._run()


class CondaPack(BasePack):
    def __init__(self, model_id, config_json):
        BasePack.__init__(self, model_id, config_json)
        self.conda = SimpleConda()
        self.logger.debug("Initializing conda packer")

    def _setup(self):
        self.logger.debug("Setting up")
        model_id = self.model_id
        installs_file = os.path.join(self._model_path(model_id), PYTHON_INSTALLS)
        self.logger.debug("Installs file {0}".format(installs_file))
        model_path = self._model_path(model_id)
        env = self.conda.specs_from_dockerfile(
            model_path, dest=None, use_checksum=USE_CHECKSUM, name=model_id
        )
        self.logger.debug("Conda environment {0}".format(env))
        if not self.conda.exists(env):
            self.logger.debug("Environment {0} does not exist".format(env))
            # clone base conda environment and add model dependencies
            base_env = self.conda.get_base_env(model_path)
            if not self.conda.exists(base_env):
                self.logger.debug(
                    "{0} base environment does not exist".format(base_env)
                )
                org = base_env.split("-")[1]
                tag = "-".join(base_env.split("-")[-2:])
                self.logger.debug("Setting up base environment {0}".format(base_env))
                SetupBaseConda().setup(org=org, tag=tag)
            self.logger.info(
                "Cloning base conda environment and adding model dependencies"
            )
            self.conda.clone(base_env, env)
            with open(installs_file, "r") as f:
                commandlines = f.read()
            self.conda.run_commandlines(environment=env, commandlines=commandlines)
        # create environment yml file
        self.logger.debug("Creating environment YAML file")
        self.conda.export_env_yml(env, model_path)
        # store conda environment in the local environment database
        self.logger.debug("Storing Conda environment in the local environment database")
        db = EnvironmentDb(config_json=self.config_json)
        db.table = "conda"
        db.insert(model_id=model_id, env=env)
        self.logger.debug("Done with the Conda setup")
        return env

    def _run(self):
        model_id = self.model_id
        env = self._setup()
        pack_snippet = """
        python {0}
        """.format(
            self.cfg.HUB.PACK_SCRIPT
        )
        self.logger.debug("Using environment {0}".format(env))
        self.logger.debug("Running command: {0}".format(pack_snippet.strip()))
        self.conda.run_commandlines(environment=env, commandlines=pack_snippet)
        self._bentoml_bundle_symlink()

    def run(self):
        self.logger.debug("Packing model with Conda")
        self._write_python_installs()
        self._run()


class DockerPack(BasePack):
    def __init__(self, model_id, config_json):
        BasePack.__init__(self, model_id, config_json)
        self.docker = SimpleDocker()
        self.docker_org = self.cfg.EXT.DOCKERHUB_ORG
        self.docker_tag = self.cfg.ENV.DOCKER.REPO_TAG
        self.logger.debug("Initializing docker packer")

    def run(self):
        self.logger.debug("Packing model with Docker")
        model_id = self.model_id
        # build docker image
        self.docker.build(".", self.docker_org, model_id, self.docker_tag)
        # pack with the docker script
        name = self.docker.run(self.docker_org, model_id, self.docker_tag, name=None)
        self.docker.exec_container(name, "python %s" % self.cfg.HUB.PACK_SCRIPT)
        # copy bundle from docker image to host
        tmp_dir = tempfile.mkdtemp()
        self.docker.cp_from_container(
            name, "/root/bentoml/repository/%s" % model_id, tmp_dir
        )
        #  save as bentoml
        mdl = bentoml_load(tmp_dir)
        mdl.save()
        self._bentoml_bundle_symlink()


def get_runner(pack_mode):
    if pack_mode == "system":
        return SystemPack
    if pack_mode == "venv":
        return VenvPack
    if pack_mode == "conda":
        return CondaPack
    if pack_mode == "docker":
        return DockerPack
