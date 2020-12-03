"""Fetch model and save it to a BentoML bundle"""

import os
import shutil
import runpy
import sys
import subprocess
import textwrap
import tempfile
from bentoml import load as bentoml_load
from .. import ErsiliaBase
from ..utils.download import GitHubDownloader, OsfDownloader, PseudoDownloader
from ..utils.zip import Zipper
from ..utils.docker import SimpleDocker
from ..utils.conda import SimpleConda
from ..utils.identifiers import FileIdentifier
from ..utils.terminal import run_command
from ..default import CONDA_ENV_YML_FILE


CONDA_INSTALLS = "conda_installs.sh"

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
        dev_path = self.cred.LOCAL.DEVEL_MODELS_PATH
        if dev_path is None:
            return None
        path = os.path.join(dev_path, model_id)
        if not os.path.exists(path):
            return None
        return path

    def _data_path(self, model_id):
        return os.path.join(self.cfg.LOCAL.DATA, model_id)

    def _get_dockerfile(self, model_id):
        fn = os.path.join(self._model_path(model_id), "Dockerfile")
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

    @staticmethod
    def _get_bentoml_location(model_id):
        cmd = ["bentoml", "get", "%s:latest" % model_id, "--print-location", "--quiet"]
        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        result = result.stdout.decode("utf-8").rstrip()
        return result

    def _get_bundle_location(self, model_id):
        return os.path.join(self._bundles_dir, model_id, self._get_latest_bundle_tag(model_id))

    def get_repo(self, model_id):
        """Fetch model from local development folder or download from GitHub"""
        folder = self._model_path(model_id)
        dev_model_path = self._dev_model_path(model_id)
        if dev_model_path is not None:
            shutil.copytree(dev_model_path, folder)
        else:
            self.github_down.clone(org=self.org, repo=model_id, destination=folder)

    def get_model(self, model_id):
        """Create a ./model folder in the model repository"""
        folder = os.path.join(self._model_path(model_id), "model")
        if os.path.exists(folder):
            return
        if self.local:
            path = os.path.join(self._data_path(model_id), "model")
            self.pseudo_down.fetch(path, folder)
        else:
            path = os.path.join("models", model_id+".zip")
            self.osf_down.fetch(project_id=self.cfg.EXT.OSF_PROJECT, filename=path,
                                destination=self._dest_dir, tmp_folder=self._tmp_dir)
            self.zipper.unzip(os.path.join(self._dest_dir, model_id+".zip"), os.path.join(self._tmp_dir, model_id))
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

    def pack(self, model_id):
        """Pack model. Greatly inspired by BentoML."""
        folder = self._model_path(model_id)
        sys.path.insert(0, folder)
        cwd = os.getcwd()
        os.chdir(folder)
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

    def as_bentoml(self, model_id):
        """Save in the system BentoML folder"""
        mdl = bentoml_load(self._get_bundle_location(model_id))
        mdl.save()

    def pip_install(self, model_id):
        """Install the model and distribute as a python package"""
        bento = self._get_bundle_location(model_id)
        run_command([sys.executable, "-m", "pip", "install", bento], quiet=True)

    def fetch(self, model_id, pip=True, bentoml=False):
        if self.overwrite:
            self.deleter.delete(model_id)
        self.get_repo(model_id)
        self.get_model(model_id)
        self.pack(model_id)
        if bentoml:
            self.as_bentoml(model_id)
        if pip:
            self.pip_install(model_id)
