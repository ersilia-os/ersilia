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
from ..utils.aux.conda_env_resolve import checksum_from_conda_yaml_file
from ..utils.terminal import run_command


class ModelFetcher(ErsiliaBase):

    def __init__(self, config_json=None, overwrite=True, local=False):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.token = self.cfg.HUB.TOKEN
        self.org = self.cfg.HUB.ORG
        self.tag = self.cfg.HUB.TAG
        self.docker_org = self.cfg.EXT.DOCKERHUB_ORG
        self.docker_tag = "repo"
        self.overwrite = overwrite
        self.pseudo_down = PseudoDownloader(overwrite=overwrite)
        self.osf_down = OsfDownloader(overwrite=overwrite)
        self.github_down = GitHubDownloader(self.token) # TODO: add overwrite?
        self.zipper = Zipper(remove=True) # TODO: Add overwrite?
        self.docker = SimpleDocker()
        self.conda = SimpleConda()
        self.local = local

    def _model_path(self, model_id):
        folder = os.path.join(self._dest_dir, model_id)
        return folder

    def _data_path(self, model_id):
        return os.path.join(self.cfg.LOCAL.DATA, model_id)

    def _get_dockerfile(self, model_id):
        fn = os.path.join(self._model_path(model_id), "Dockerfile")
        if os.path.exists(fn):
            return fn
        else:
            return None

    def _get_conda_env_yaml_file(self, model_id):
        fn = os.path.join(self._model_path(model_id), "environment.yml")
        if os.path.exists(fn):
            return fn
        else:
            return None

    @staticmethod
    def _get_bentoml_location(model_id):
        cmd = ["bentoml", "get", "%s:latest" % model_id, "--print-location", "--quiet"]
        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        result = result.stdout.decode("utf-8").rstrip()
        return result

    def _get_bundle_location(self, model_id):
        return os.path.join(self._bundles_dir, model_id, self._get_latest_bundle_tag(model_id))

    def get_repo(self, model_id):
        """Download the model from GitHub"""
        folder = self._model_path(model_id)
        self.github_down.clone(org=self.org, repo=model_id, destination=folder)

    def get_model(self, model_id):
        """Create a ./model folder in the model repository"""
        folder = os.path.join(self._model_path(model_id), "model")
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

    def _run_pack_script(self, folder):
        script_path = os.path.join(folder, self.cfg.HUB.PACK_SCRIPT)
        runpy.run_path(script_path)

    def _run_pack_script_conda(self, folder):
        yml = self._get_conda_env_yaml_file(folder)
        checksum = checksum_from_conda_yaml_file(yml, overwrite=True)
        default_env = self.conda.default_env()
        is_base = self.conda.is_base()
        if not is_base:
            bash_script = """
            source ${0}/etc/profile.d/conda.sh
            conda deactivate
            """.format(self.conda.conda_prefix(False))
        else:
            bash_script = ""
        bash_script += """
        source ${0}/etc/profile.d/conda.sh
        cd {1}
        """.format(self.conda.conda_prefix(True), folder)
        if not self.conda.exists(checksum):
            bash_script += """
            conda env create -f environment.yml
            """
        bash_script += """
        conda activate {0}
        python {1}
        conda deactivate""".format(checksum, self.cfg.HUB.PACK_SCRIPT)
        if not self.conda.is_base():
            bash_script += """
            conda activate {0}
            """.format(default_env)
        bash_script = textwrap.dedent(bash_script)
        tmp_folder = tempfile.mkdtemp()
        fn = os.path.join(tmp_folder, "create_env.sh")
        with open(fn, "w") as f:
            f.write(bash_script)
        with open(fn, "r") as f:
            print(f.read())
        run_command("bash {0}".format(fn), quiet=True)

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
        """Pack model"""
        folder = self._model_path(model_id)
        sys.path.insert(0, folder)
        cwd = os.getcwd()
        os.chdir(folder)
        yf = self._get_conda_env_yaml_file(model_id)
        df = self._get_dockerfile(model_id)
        if yf is not None:
            # pack using conda
            self._run_pack_script_conda(folder)
            sys.exit()
        elif df is not None:
            # pack using docker
            self._run_pack_script_docker(model_id)
        else:
            # try to run without conda environment or docker (not recommended)
            self._run_pack_script(folder)
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
        self.get_repo(model_id)
        self.get_model(model_id)
        self.pack(model_id)
        if bentoml:
            self.as_bentoml(model_id)
        if pip:
            self.pip_install(model_id)
