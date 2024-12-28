import os

from packaging import version

from .. import logger
from ..utils.conda import SimpleConda
from ..utils.logging import make_temp_dir
from ..utils.terminal import run_command
from ..utils.versioning import Versioner
from .utils.clone import ErsiliaCloner


class SetupBaseConda(object):
    """
    A class to handle the setup of a base Conda environment for the model server.

    Methods
    -------
    find_closest_python_version(python_version)
        Finds the closest available Python version.
    setup(org, tag)
        Creates a Conda environment to be used as base environment for the model server.
    delete(org, tag)
        Deletes the specified Conda environment.
    """

    def __init__(self, config_json=None):
        self.conda = SimpleConda()
        self.versions = Versioner()
        self.cloner = ErsiliaCloner(config_json=config_json)

    @staticmethod
    def _is_ersiliaos(org: str) -> bool:
        """
        Checks if the organization is ErsiliaOS.

        Parameters
        ----------
        org : str
            The organization name.

        Returns
        -------
        bool
            True if the organization is ErsiliaOS, False otherwise.
        """
        if org == "ersiliaos":
            return True
        else:
            return False

    @staticmethod
    def _is_bentoml(org: str) -> bool:
        """
        Checks if the organization is BentoML.

        Parameters
        ----------
        org : str
            The organization name.

        Returns
        -------
        bool
            True if the organization is BentoML, False otherwise.
        """
        if org == "bentoml":
            return True
        else:
            return False

    def _parse_tag(self, tag: str) -> dict:
        """
        Parses the tag to extract version and Python version information.

        Parameters
        ----------
        tag : str
            The tag in the format '0.0.0-py37'.

        Returns
        -------
        dict
            A dictionary containing version and Python version information.
        """
        tag = tag.split("-")
        data = {
            "ver": tag[0],
            "py": tag[1],
            "python": self.versions.reformat_py(tag[1]),
        }
        return data

    def _install_command(self, org: str, tag: str) -> str:
        """
        Generates the install command based on the organization and tag.

        Parameters
        ----------
        org : str
            The organization name.
        tag : str
            The tag in the format '0.0.0-py37'.

        Returns
        -------
        str
            The install command.
        """
        tag = self._parse_tag(tag)
        if self._is_bentoml(org):
            if tag["ver"] == "0.11.0":
                logger.debug("Installing from ersilia's custom BentoML")
                cmd = "python -m pip install git+https://github.com/ersilia-os/bentoml-ersilia.git"
            else:
                logger.debug("Installing from BentoML directly")
                cmd = "python -m pip install bentoml=={0}".format(tag["ver"])
        elif self._is_ersiliaos(org):
            cmd = "python -m pip install -e ."
        else:
            raise Exception
        return cmd

    def _get_env_name(self, org: str, tag: str) -> str:
        """
        Generates the environment name based on the organization and tag.

        Parameters
        ----------
        org : str
            The organization name.
        tag : str
            The tag in the format '0.0.0-py37'.

        Returns
        -------
        str
            The environment name.
        """
        return self.versions.base_conda_name(org, tag)

    def find_closest_python_version(self, python_version: str) -> str:
        """
        Finds the closest available Python version.

        Parameters
        ----------
        python_version : str
            The desired Python version.

        Returns
        -------
        str
            The closest available Python version.
        """
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "conda_search_python.txt")
        tmp_script = os.path.join(tmp_folder, "script.sh")
        is_base = self.conda.is_base()
        bash_script = """
        source {0}/etc/profile.d/conda.sh
        conda search python > {1}
        """.format(self.conda.conda_prefix(is_base), tmp_file)
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        run_command("bash {0}".format(tmp_script))
        with open(tmp_file, "r") as f:
            available_versions = []
            for l in f:
                if l.startswith("python"):
                    available_versions += [
                        l.split("python")[1].lstrip(" ").split(" ")[0]
                    ]
        kept = []
        for v in available_versions:
            if version.parse(v) >= version.parse(python_version):
                kept += [v]
        return ".".join(kept[0].split(".")[:2])

    def setup(self, org: str, tag: str) -> None:
        """
        Creates a Conda environment to be used as base environment for the model server.

        Parameters
        ----------
        org : str
            The organization name (bentoml or ersiliaos).
        tag : str
            The tag in the format '0.0.0-py37'.

        Returns
        -------
        None
        """
        env = self._get_env_name(org, tag)
        if self.conda.exists(env):
            return
        ptag = self._parse_tag(tag)
        cmd = self._install_command(org, tag)
        tmp_folder = make_temp_dir(prefix="ersilia-")
        if self._is_ersiliaos(org):
            tmp_repo = self.cloner.clone(tmp_folder, version=ptag["ver"])
        else:
            tmp_repo = tmp_folder
        tmp_script = os.path.join(tmp_folder, "script.sh")
        is_base = self.conda.is_base()
        if not is_base:
            bash_script = """
            source {0}/etc/profile.d/conda.sh
            conda deactivate
            """.format(self.conda.conda_prefix(False))
        else:
            bash_script = ""
        bash_script += """
        source {0}/etc/profile.d/conda.sh
        """.format(self.conda.conda_prefix(True))
        python_version = self.find_closest_python_version(ptag["python"])
        bash_script += """
        cd {0}
        conda create --no-default-packages -n {1} python={2} -y
        conda activate {1}
        {3}
        conda deactivate
        """.format(tmp_repo, env, python_version, cmd)
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        run_command("bash {0}".format(tmp_script))

    def delete(self, org: str, tag: str) -> None:
        """
        Deletes the specified Conda environment.

        Parameters
        ----------
        org : str
            The organization name.
        tag : str
            The tag in the format '0.0.0-py37'.

        Returns
        -------
        None
        """
        env = self._get_env_name(org, tag)
        self.conda.delete(env)
