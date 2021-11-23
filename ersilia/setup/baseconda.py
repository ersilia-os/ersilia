import tempfile
import sys
import os

from ..utils.conda import SimpleConda
from ..utils.terminal import run_command
from ..utils.versioning import Versioner
from .utils.clone import ErsiliaCloner


class SetupBaseConda(object):
    def __init__(self, config_json=None):
        self.conda = SimpleConda()
        self.versions = Versioner()
        self.cloner = ErsiliaCloner(config_json=config_json)

    @staticmethod
    def _is_ersiliaos(org):
        if org == "ersiliaos":
            return True
        else:
            return False

    @staticmethod
    def _is_bentoml(org):
        if org == "bentoml":
            return True
        else:
            return False

    def _parse_tag(self, tag):
        tag = tag.split("-")
        return {
            "ver": tag[0],
            "py": tag[1],
            "python": self.versions.reformat_py(tag[1]),
        }

    def _install_command(self, org, tag):
        tag = self._parse_tag(tag)
        if self._is_bentoml(org):
            cmd = "pip install bentoml=={0}".format(tag["ver"])
        elif self._is_ersiliaos(org):
            cmd = "pip install -e ."
        else:
            raise Exception
        return cmd

    def _get_env_name(self, org, tag):
        return self.versions.base_conda_name(org, tag)

    def setup(self, org, tag):
        """Creates a conda enviornment to be used as base environment to be used as model server.

        Args:
            org: organisation (bentoml or ersiliaos)
            tag: 0.0.0-py37 format
        """
        env = self._get_env_name(org, tag)
        if self.conda.exists(env):
            return
        ptag = self._parse_tag(tag)
        cmd = self._install_command(org, tag)
        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        if self._is_ersiliaos(org):
            tmp_repo = self.cloner.clone(tmp_folder, version=ptag["ver"])
        else:
            tmp_repo = tmp_folder
        tmp_script = os.path.join(tmp_folder, "script.sh")
        is_base = self.conda.is_base()
        if not is_base:
            bash_script = """
            source ${0}/etc/profile.d/conda.sh
            conda deactivate
            """.format(
                self.conda.conda_prefix(False)
            )
        else:
            bash_script = ""
        bash_script += """
        source ${0}/etc/profile.d/conda.sh
        """.format(
            self.conda.conda_prefix(True)
        )
        bash_script += """
        cd {0}
        conda create --no-default-packages -n {1} python={2} -y
        conda activate {1}
        {3}
        conda deactivate
        """.format(
            tmp_repo, env, ptag["python"], cmd
        )
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        run_command("bash {0}".format(tmp_script))

    def delete(self, org, tag):
        env = self._get_env_name(org, tag)
        self.conda.delete(env)
