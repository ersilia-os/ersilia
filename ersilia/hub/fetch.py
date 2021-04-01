"""Fetch model and save it to a BentoML bundle"""

import os
import shutil
import runpy
import sys
import textwrap
import tempfile
import importlib
import yaml
from bentoml import load as bentoml_load
from .. import ErsiliaBase
from ..utils.download import GitHubDownloader, OsfDownloader, PseudoDownloader
from ..utils.zip import Zipper
from ..utils.docker import SimpleDocker
from ..utils.conda import SimpleConda
from ..utils.identifiers.file import FileIdentifier
from ..utils.terminal import run_command
from ..utils.paths import Paths
from ..default import CONDA_ENV_YML_FILE
from ..db.environments.localdb import EnvironmentDb
from .repo import ServiceFile, DockerfileFile


CONDA_INSTALLS = "conda_installs.sh"

DOCKERFILE = "Dockerfile"
ENVIRONMENT_YML = "environment.yml"


class ModelStatus(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def is_downloaded(self, model_id):
        essentials = ["README.md", "model"] # essential files
        dst_dir = os.path.join(self._dest_dir, model_id)
        if not os.path.exists(dst_dir):
            return False
        items = {i for i in os.listdir(dst_dir)}
        for essential in essentials:
            if essential not in items:
                return False
        return True

    def is_docker(self, model_id):
        db = EnvironmentDb(config_json=self.config_json)
        db.table = "docker"
        docker = SimpleDocker()
        envs = db.envs_of_model(model_id)
        for env in envs:
            img, tag = env.split(":")
            if docker.exists(self.cfg.EXT.DOCKERHUB_ORG, img, tag):
                return True
        return False

    def is_conda(self, model_id):
        db = EnvironmentDb(config_json=self.config_json)
        db.table = "conda"
        conda = SimpleConda()
        envs = db.envs_of_model(model_id)
        for env in envs:
            if conda.exists(env):
                return True
        return False

    def is_pip(self, model_id):
        try:
            importlib.import_module(model_id, package=None)
            return True
        except ModuleNotFoundError:
            return False

    def _is_bento_folder(self, model_folder):
        if not os.path.exists(model_folder):
            return False
        essentials = ["bentoml.yml"]
        items = {i for i in os.listdir(model_folder)}
        for essential in essentials:
            if essential not in items:
                return False
        return True

    def is_bundle(self, model_id):
        model_folder = self._get_bundle_location(model_id)
        return self._is_bento_folder(model_folder)

    def is_bentoml(self, model_id):
        model_folder = self._get_bentoml_location(model_id)
        return self._is_bento_folder(model_folder)

    def status(self, model_id):
        """Check installation and deployment status of the model"""
        results = {
            "download": self.is_downloaded(model_id),
            "bentoml": self.is_bentoml(model_id),
            "bundle": self.is_bundle(model_id),
            "docker": self.is_docker(model_id),
            "conda": self.is_conda(model_id),
            "pip": self.is_pip(model_id)
        }
        return results


class ModelFetcher(ErsiliaBase):

    def __init__(self,
                 config_json=None, credentials_json=None,
                 overwrite=True, local=False):
        ErsiliaBase.__init__(self,
                             config_json=config_json,
                             credentials_json=credentials_json)
        self.token = self.cfg.HUB.TOKEN
        self.org = self.cfg.HUB.ORG
        self.tag = self.cfg.HUB.TAG
        self.docker_org = self.cfg.EXT.DOCKERHUB_ORG
        self.docker_tag = self.cfg.ENV.DOCKER.REPO_TAG
        self.overwrite = overwrite
        self.pseudo_down = PseudoDownloader(overwrite=overwrite)
        self.osf_down = OsfDownloader(overwrite=overwrite)
        self.github_down = GitHubDownloader(self.token) # TODO: add overwrite?
        self.zipper = Zipper(remove=True) # TODO: Add overwrite?
        self.docker = SimpleDocker()
        self.conda = SimpleConda()
        self.file_identifier = FileIdentifier()
        self.local = local
        if self.overwrite:
            from .delete import ModelFullDeleter
            self.deleter = ModelFullDeleter(config_json=config_json)
        else:
            self.deleter = None

    def _file_checksum(self, path):
        return self.file_identifier.encode(path, n=CHECKSUM_NCHAR)

    def _model_path(self, model_id):
        folder = os.path.join(self._dest_dir, model_id)
        return folder

    def _dev_model_path(self, model_id):
        if not self.cred.exists:
            return None
        dev_path = self.cred.LOCAL.DEVEL_MODELS_PATH
        if dev_path is None:
            return None
        path = os.path.join(dev_path, model_id)
        if os.path.exists(path):
            return path
        else:
            path = Paths().ersilia_development_path()
            path = os.path.join(path, "test", "models", model_id)
            if os.path.exists(path):
                return path
        return None

    def _data_path(self, model_id):
        return os.path.join(self.cfg.LOCAL.DATA, model_id)

    def _get_dockerfile(self, model_id):
        fn = os.path.join(self._model_path(model_id), DOCKERFILE)
        if os.path.exists(fn):
            return fn
        else:
            return None

    def _get_conda_env_yaml_file(self, model_id):
        fn = os.path.join(self._model_path(model_id), CONDA_ENV_YML_FILE)
        if os.path.exists(fn):
            return fn
        else:
            return None

    def _get_conda_installs_from_dockerfile(self, model_id):
        fn = self._get_dockerfile(model_id)
        if fn is None:
            return None
        runs = self.conda.get_install_commands_from_dockerfile(fn)
        if not runs:
            return None
        fn = os.path.join(self._model_path(model_id), CONDA_INSTALLS)
        with open(fn, "w") as f:
            for r in runs:
                f.write("{0}\n".format(r))
        return fn

    def get_repo(self, model_id):
        """Fetch model from local development folder or download from GitHub"""
        folder = self._model_path(model_id)
        dev_model_path = self._dev_model_path(model_id)
        if dev_model_path is not None:
            shutil.copytree(dev_model_path, folder)
        else:
            self.github_down.clone(org=self.org, repo=model_id, destination=folder)

    def get_model_parameters(self, model_id):
        """Create a ./model folder in the model repository"""
        folder = os.path.join(self._model_path(model_id), "model")
        if os.path.exists(folder):
            return
        if self.local:
            path = os.path.join(self._data_path(model_id), "model")
            self.pseudo_down.fetch(path, folder)
        else:
            path = os.path.join("models", model_id+".zip")
            try:
                self.osf_down.fetch(project_id=self.cfg.EXT.OSF_PROJECT, filename=path,
                                    destination=self._dest_dir, tmp_folder=self._tmp_dir)
                zip_file = os.path.join(self._dest_dir, model_id+".zip")
                self.zipper.unzip(zip_file, os.path.join(self._tmp_dir, model_id))
                src = os.path.join(self._tmp_dir, model_id)
                dst = os.path.join(self._dest_dir, model_id, "model")
                if os.path.exists(dst):
                    if self.overwrite:
                        shutil.rmtree(dst)
                    else:
                        shutil.rmtree(src)
                        return
                shutil.move(os.path.join(src, "model"), os.path.join(self._dest_dir, model_id))
                shutil.rmtree(src)
            except:
                if not os.path.exists(folder):
                    os.mkdir(folder)

    def _repath_pack_save(self, model_id):
        folder = self._model_path(model_id)
        path = os.path.join(folder, self.cfg.HUB.PACK_SCRIPT)
        with open(path, "r") as f:
            text = f.read()
        text = text.replace(
            "service.save()",
            "service.save('{0}')".format(self._bundles_dir)
        )
        with open(path, "w") as f:
            f.write(text)

    def _run_pack_script(self, model_id):
        self._repath_pack_save(model_id)
        folder = self._model_path(model_id)
        script_path = os.path.join(folder, self.cfg.HUB.PACK_SCRIPT)
        runpy.run_path(script_path)

    def _setup_conda(self, model_id):
        installs_file = os.path.join(self._model_path(model_id), CONDA_INSTALLS)
        model_path = self._model_path(model_id)
        checksum = self.conda.checksum_from_dockerfile(model_path)
        if not self.conda.exists(checksum):
            # clone base conda environment and add model dependencies
            base_env = self.conda.get_base_env(model_path)
            self.conda.clone(base_env, checksum)
            with open(installs_file, "r") as f:
                commandlines = f.read()
            self.conda.run_commandlines(environment=checksum, commandlines=commandlines)
        # create environment yml file
        self.conda.export_env_yml(checksum, model_path)
        # store conda environment in the local environment database
        db = EnvironmentDb(config_json=self.config_json)
        db.table = "conda"
        db.insert(model_id=model_id, env=checksum)
        return checksum

    def _run_pack_script_conda(self, model_id):
        checksum = self._setup_conda(model_id)
        self._repath_pack_save(model_id)
        pack_snippet = """
        python {0}
        """.format(self.cfg.HUB.PACK_SCRIPT)
        self.conda.run_commandlines(environment=checksum, commandlines=pack_snippet)

    def _run_pack_script_docker(self, model_id):
        # build docker image
        self.docker.build(".", self.docker_org, model_id, self.docker_tag)
        # pack with the docker script
        name = self.docker.run(self.docker_org, model_id, self.docker_tag, name=None)
        self.docker.exec_container(name, "python %s" % self.cfg.HUB.PACK_SCRIPT)
        # copy bundle from docker image to host
        self.docker.cp_from_container(name,
                                      "/root/bentoml/repository/%s" % model_id,
                                      self._bundles_dir)

    def _bundle_uses_ersilia(self, model_id):
        """Check if the bundle imports ersilia"""
        src = os.path.join(self._get_bundle_location(model_id), model_id, "src")
        tmp_folder = tempfile.mkdtemp()
        tmp_file = os.path.join(tmp_folder, "grep.txt")
        cmd = "grep -R 'ersilia' */* > {0}".format(tmp_file)
        run_command(cmd, quiet=True)
        with open(tmp_file, "r") as f:
            grep = f.read()
        if grep:
            return True
        else:
            return False

    def _modify_bundle_environment_yml(self, model_id):
        """This method corrects some inconsistencies in default yml file created by conda.
        - libfortran version is removed (inconsistency between mac and linux)
        - ruamel is not explicitly installed
        - ersilia is downloaded from github
        """
        dir = self._get_bundle_location(model_id)
        yml_file = os.path.join(dir, ENVIRONMENT_YML)
        if not os.path.exists(yml_file):
            return
        try:
            with open(yml_file, "r") as f:
                data = yaml.safe_load(f)
                for i, d in enumerate(data["dependencies"][:-1]):
                    if "libgfortran=" in d:
                        data["dependencies"][i] = "libgfortran"
                for i, p in enumerate(data["dependencies"][-1]["pip"]):
                    if "ruamel" in p:
                        data["dependencies"][-1]["pip"][i] = None
                v = [x for x in data["dependencies"][-1]["pip"] if x is not None]
                data["dependencies"][-1]["pip"] = v
            with open(yml_file, "w") as f:
                yaml.safe_dump(data, f)
        except:
            return

    def _bundle_environment_yml_has_ersilia(self, model_id):
        """Check if bundle environment.yml file installs ersilia"""
        search = "ersilia="
        dir = self._get_bundle_location(model_id)
        yml_file = os.path.join(dir, ENVIRONMENT_YML)
        if not os.path.exists(yml_file):
            return None
        with open(yml_file, "r") as f:
            data = yaml.safe_load(f)
            for i, d in enumerate(data["dependencies"][:-1]):
                if search in d:
                    return True
            for i, p  in enumerate(data["dependencies"][-1]["pip"]):
                if search in p:
                    return True
        return False

    def _bundle_dockerfile_has_ersilia(self, model_id):
        """Check if bundle Dockerfile uses ersilia"""
        dir = self._get_bundle_location(model_id)
        dockerfile = os.path.join(dir, DOCKERFILE)
        if not os.path.exists(dockerfile):
            return None
        tmp_folder = tempfile.mkdtemp()
        tmp_file = os.path.join(tmp_folder, "grep.txt")
        cmd = "grep -R 'ersilia' {0} > {1}".format(dockerfile, tmp_file)
        run_command(cmd, quiet=True)
        with open(tmp_file, "r") as f:
            grep = f.read()
        if grep:
            return True
        else:
            return False

    def _modify_bundle_dockerfile(self, model_id):
        """This method modifies the dockerfile generated by BentoML.
        If ersilia is needed by the bundle and not specified, it is added in the Dockerfile.
        If a development path exists locally, then this is the one we use.
        If not, we use the latest version available from github. This may change in the future to use the PyPi version, correspondingly.
        """

        if not self._bundle_uses_ersilia(model_id):
            return
        if self._bundle_environment_yml_has_ersilia(model_id):
            return
        if self._bundle_dockerfile_has_ersilia(model_id):
            return
        dockerfile = os.path.join(self._get_bundle_location(model_id), DOCKERFILE)
        text = ["", "# Install ersilia"]
        text += ["RUN pip install git+https://github.com/{0}/{1}.git".format(self.cfg.HUB.ORG, self.cfg.HUB.PACKAGE)] # TO-DO add version with the @ character
        with open(dockerfile, "r") as f:
            lines = []
            for l in f:
                lines += [l.rstrip()]
        if lines[1][:10] == "MAINTAINER":
            lines = lines[:2] + text + lines[2:]
        else:
            lines = lines[:1] + text + lines[1:]
        with open(dockerfile, "w") as f:
            for l in lines:
                f.write(l + os.linesep)

    def pack(self, model_id):
        """Pack model. Greatly inspired by BentoML."""
        folder = self._model_path(model_id)
        ServiceFile(folder).rename_service()
        sys.path.insert(0, folder)
        cwd = os.getcwd()
        os.chdir(folder)
        dockerfile = DockerfileFile(folder)
        # if dockerfile has no runs (i.e. no installation is necessary, apart from bentoml),
        #  run pack script from the current environment
        if not dockerfile.has_runs():
            # get bentoml version
            res = dockerfile.get_bentoml_version()
            if res is None:
                raise Exception
            bentoml_version = res["version"]
            # TODO: use python virtual environment to use exact bentoml version
            self._run_pack_script(model_id)
        else:
            cf = self._get_conda_installs_from_dockerfile(model_id)
            df = self._get_dockerfile(model_id)
            if cf is not None:
                # pack using conda
                self._run_pack_script_conda(model_id)
            elif df is not None:
                # pack using docker
                self._run_pack_script_docker(model_id)
            else:
                # try to run without conda environment or docker (not recommended)
                self._run_pack_script(model_id)
        os.chdir(cwd)
        sys.path.remove(folder)
        # Slightly modify bundle environment YAML file, if exists
        self._modify_bundle_environment_yml(model_id)
        # Slightly modify bundle Dockerfile, if necessary.
        self._modify_bundle_dockerfile(model_id)

    def as_bentoml(self, model_id):
        """Save in the system BentoML folder"""
        mdl = bentoml_load(self._get_bundle_location(model_id))
        mdl.save()

    def pip_install(self, model_id):
        """Install the model and distribute as a python package"""
        bento = self._get_bundle_location(model_id)
        run_command([sys.executable, "-m", "pip", "install", bento], quiet=True)

    def containerize(self, model_id):
        """Containerize model using bentoml"""
        bento = self._get_bundle_location(model_id)
        tag = self.cfg.ENV.DOCKER.LATEST_TAG
        bundle_tag = self._get_latest_bundle_tag(model_id)
        run_command("bentoml containerize {1}:{2} -t {0}/{1}:{2}".format(self.docker_org, model_id, tag), quiet=True)
        # store docker in the local environment database
        db = EnvironmentDb(config_json=self.config_json)
        db.table = "docker"
        db.insert(model_id=model_id, env="{0}/{1}:{2}".format(self.docker_org, model_id, tag))

    def fetch(self, model_id, pip=True, containerize=False, bentoml=False):
        if self.overwrite:
            self.deleter.delete(model_id)
        ms = ModelStatus()
        if not ms.is_downloaded(model_id):
            self.get_repo(model_id)
            self.get_model_parameters(model_id)
            self.pack(model_id)
        if bentoml:
            if not ms.is_bentoml(model_id):
                self.as_bentoml(model_id)
        if containerize:
            if not ms.is_docker(model_id):
                self.containerize(model_id)
        if pip:
            if not ms.is_pip(model_id):
                self.pip_install(model_id)
