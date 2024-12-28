import os

try:
    import bentoml
except:
    bentoml = None

from ..... import throw_ersilia_exception
from .....db.environments.localdb import EnvironmentDb
from .....db.environments.managers import DockerManager
from .....default import DEFAULT_VENV
from .....setup.baseconda import SetupBaseConda
from .....utils.conda import SimpleConda
from .....utils.docker import SimpleDocker
from .....utils.exceptions_utils.fetch_exceptions import CondaEnvironmentExistsError
from .....utils.logging import make_temp_dir
from .....utils.terminal import run_command
from .....utils.venv import SimpleVenv
from ... import MODEL_INSTALL_COMMANDS_FILE
from . import BasePack

USE_CHECKSUM = False


class SystemPack(BasePack):
    """
    A class used to pack models with system installation.

    Methods
    -------
    run()
        Run the system packer.

    Examples
    --------
    .. code-block:: python

        packer = SystemPack(
            model_id="eosxxxx", config_json=config
        )
        packer.run()
    """

    def __init__(self, model_id, config_json):
        BasePack.__init__(self, model_id, config_json)
        self.logger.debug("Initializing system packer")

    def _run(self):
        self.logger.debug("Packing model with system installation")
        folder = self._model_path(self.model_id)
        script_path = os.path.join(folder, self.cfg.HUB.PACK_SCRIPT)
        run_command("python {0}".format(script_path))
        self._symlinks()

    def run(self):
        """
        Run the system packer.
        """
        self._run()


class VenvPack(BasePack):
    """
    A class used to pack models with virtual environments.

    Methods
    -------
    run()
        Run the virtual environment packer.
    """

    def __init__(self, model_id, config_json):
        BasePack.__init__(self, model_id, config_json)
        self.logger.debug("Initializing virtualenv packer")

    def _setup(self) -> SimpleVenv:
        model_path = self._model_path(self.model_id)
        installs_file = os.path.join(model_path, MODEL_INSTALL_COMMANDS_FILE)
        self.logger.debug("Reading install commands from {0}".format(installs_file))
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
        """.format(self.cfg.HUB.PACK_SCRIPT)
        venv.run_commandlines(environment=DEFAULT_VENV, commandlines=pack_snippet)
        self._symlinks()

    def run(self):
        """
        Run the virtual environment packer.
        """
        self.logger.debug("Packing model with VirtualEnv")
        self._write_model_install_commands()
        self._run()


class CondaPack(BasePack):
    """
    A class used to pack models with Conda environments.

    Methods
    -------
    run()
        Run the Conda environment packer.
    """

    def __init__(self, model_id, config_json):
        BasePack.__init__(self, model_id, config_json)
        self.conda = SimpleConda()
        self.logger.debug("Initializing conda packer")

    def _setup(self) -> str:
        self.logger.debug("Setting up")
        model_id = self.model_id
        installs_file = os.path.join(
            self._model_path(model_id), MODEL_INSTALL_COMMANDS_FILE
        )
        self.logger.debug("Installs file {0}".format(installs_file))
        model_path = self._model_path(model_id)
        env = self.conda.specs_from_dockerfile(
            model_path, dest=None, use_checksum=USE_CHECKSUM, name=model_id
        )
        self.logger.debug("Conda environment {0}".format(env))
        if not self.conda.exists(env):
            self.logger.debug("Environment {0} does not exist".format(env))
            base_env = self.conda.get_base_env(model_path)
            if "-slim-" in base_env:
                self.logger.debug("Removing -slim- from base environment!")
                base_env = base_env.replace("-slim", "")
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
            if base_env is not None:
                if not self.conda.exists(base_env):
                    raise CondaEnvironmentExistsError(base_env)

            self.conda.clone(base_env, env)
            with open(installs_file, "r") as f:
                commandlines = f.read()
            self.conda.run_commandlines(environment=env, commandlines=commandlines)
        else:
            self.logger.debug("Environment {0} does exist".format(env))
        self.logger.debug("Creating environment YAML file")
        self.conda.export_env_yml(env, model_path)
        self.logger.debug("Storing Conda environment in the local environment database")
        db = EnvironmentDb(config_json=self.config_json)
        db.table = "conda"
        db.insert(model_id=model_id, env=env)
        self.logger.debug("Done with the Conda setup")

        if env is not None:
            if not self.conda.exists(env):
                raise CondaEnvironmentExistsError(env)
        return env

    def _run(self):
        env = self._setup()
        pack_snippet = """
        python {0}
        """.format(self.cfg.HUB.PACK_SCRIPT)
        self.logger.debug("Using environment {0}".format(env))
        self.logger.debug("Running command: {0}".format(pack_snippet.strip()))
        self.conda.run_commandlines(environment=env, commandlines=pack_snippet)
        self.logger.debug(
            "Previous command successfully run inside {0} conda environment".format(env)
        )
        self.logger.debug("Now trying to establish symlinks")
        self._symlinks()

    @throw_ersilia_exception()
    def run(self):
        """
        Run the Conda environment packer.
        """
        self.logger.debug("Packing model with Conda")
        self._write_model_install_commands()
        self._run()


class DockerPack(BasePack):
    """
    A class used to pack models with Docker.

    Methods
    -------
    run()
        Run the Docker packer.
    """

    def __init__(self, model_id, config_json):
        BasePack.__init__(self, model_id, config_json)
        self.docker = SimpleDocker()
        self.docker_org = self.cfg.EXT.DOCKERHUB_ORG
        self.docker_tag = self.cfg.ENV.DOCKER.REPO_TAG
        self.logger.debug("Initializing docker packer")

    def _load_model_from_tmp(self, path: str):
        path = os.path.join(path, self.model_id)
        if not os.path.exists(path):
            return None
        items = sorted(os.listdir(path))
        if not items:
            return None
        else:
            tag = items[-1]
        path = os.path.join(path, tag)
        return bentoml.load(path)

    def _setup(self):
        name = self.docker._image_name(self.docker_org, self.model_id, self.docker_tag)
        self.logger.debug(
            "Storing Docker image {0} in the local environment database".format(name)
        )
        db = EnvironmentDb(config_json=self.config_json)
        db.table = "docker"
        db.insert(model_id=self.model_id, env=name)
        self.logger.debug("Done with the Docker setup")

    def _delete_packer_container(self):
        dm = DockerManager(config_json=self.config_json)
        dm.delete_images(self.model_id)

    def _run(self):
        model_id = self.model_id
        self._setup()
        self.logger.debug("Building docker image")
        self.docker.build(".", self.docker_org, model_id, self.docker_tag)
        self.logger.debug("Running docker")
        name = self.docker.run(self.docker_org, model_id, self.docker_tag, name=None)
        self.logger.debug("Executing container {0}".format(name))
        self.docker.exec_container(name, "python %s" % self.cfg.HUB.PACK_SCRIPT)
        self.logger.debug("Copying bundle from docker image to host")
        tmp_dir = make_temp_dir(prefix="ersilia-")
        self.logger.debug("Using this temporary directory: {0}".format(tmp_dir))
        self.docker.cp_from_container(
            name, "/root/bentoml/repository/%s" % model_id, tmp_dir
        )
        self.logger.debug("Loading bentoml")
        mdl = self._load_model_from_tmp(tmp_dir)
        mdl.save()
        self._symlinks()
        self._delete_packer_container()

    def run(self):
        """
        Run the Docker packer.
        """
        self.logger.debug("Packing model with Docker")
        self._write_model_install_commands()
        self._run()


def get_runner(pack_mode: str):
    """
    Get the runner class based on the pack mode.

    Parameters
    ----------
    pack_mode : str
        Packaging mode.

    Returns
    -------
    type
        Runner class corresponding to the pack mode.
    """
    if pack_mode == "system":
        return SystemPack
    if pack_mode == "venv":
        return VenvPack
    if pack_mode == "conda":
        return CondaPack
    if pack_mode == "docker":
        return DockerPack
