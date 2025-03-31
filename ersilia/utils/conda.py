# TODO Implement conda python telemetry for resources during model execution
import hashlib
import json
import os
import shutil
from collections import OrderedDict, defaultdict

from .. import logger, throw_ersilia_exception
from ..default import CONDA_ENV_YML_FILE
from ..utils.exceptions_utils.fetch_exceptions import ModelPackageInstallError
from ..utils.logging import make_temp_dir
from .docker import SimpleDockerfileParser
from .supp.conda_env_resolve import CHECKSUM_FILE, CHECKSUM_NCHAR
from .terminal import run_command, run_command_check_output
from .versioning import Versioner

BASE = "base"
SPECS_JSON = ".specs.json"


class BaseConda(object):
    """
    Base class for managing conda environments.

    Methods
    -------
    default_env()
        Get the default conda environment.
    is_base()
        Check if the current environment is the base environment.
    conda_prefix(is_base)
        Get the conda prefix path.
    """

    def __init__(self):
        self.SPECS_JSON = SPECS_JSON
        self.CHECKSUM_FILE = CHECKSUM_FILE

    @staticmethod
    def default_env():
        """
        Get the default conda environment.

        Returns
        -------
        str
            The name of the default conda environment.
        """
        if "CONDA_DEFAULT_ENV" in os.environ:
            return os.environ["CONDA_DEFAULT_ENV"]
        else:
            return BASE

    def is_base(self):
        """
        Check if the current environment is the base environment.

        Returns
        -------
        bool
            True if the current environment is the base environment, False otherwise.
        """
        default_env = self.default_env()
        if default_env == BASE:
            return True
        else:
            return False

    @staticmethod
    def conda_prefix(is_base):
        """
        Get the conda prefix path.

        Parameters
        ----------
        is_base : bool
            Whether the current environment is the base environment.

        Returns
        -------
        str
            The conda prefix path.
        """
        o = run_command_check_output("which conda").rstrip()
        if o:
            o = os.path.abspath(os.path.join(o, "..", ".."))
            return o
        if is_base:
            o = run_command_check_output("echo $CONDA_PREFIX").rstrip()
            return o
        else:
            o = run_command_check_output("echo $CONDA_PREFIX_1").rstrip()
            return o


class CondaUtils(BaseConda):
    """
    Utility class for managing conda environments and packages.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings in JSON format. Default is None.

    Examples
    --------
    .. code-block:: python

        conda_utils = CondaUtils()
        conda_utils.create("myenv", "3.8")

    """

    def __init__(self, config_json=None):
        BaseConda.__init__(self)
        self.versions = Versioner(config_json=config_json)

    @staticmethod
    def _parse_install(r):
        r = r.strip().split()
        tool = None
        channel = None
        packages = None
        if r[0] == "conda" and r[1] == "install" and r[2] == "-c":
            tool = "conda"
            channel = r[3]
            packages = r[4:]
        if r[0] == "conda" and r[1] == "install" and r[2] != "-c":
            tool = "conda"
            channel = "default"
            packages = r[2:]
        if (r[0] == "pip" or r[0] == "pip3") and r[1] == "install":
            tool = "pip"
            packages = r[2:]
        if tool is None:
            return None
        packages = sorted(packages)
        result = {"tool": tool, "channel": channel, "packages": packages}
        return result

    @staticmethod
    def _text_checksum(text):
        checksum = hashlib.md5(text.encode("utf-8")).hexdigest()[:CHECKSUM_NCHAR]
        return checksum

    def checksum_from_file(self, filename):
        """
        Generate a checksum for the contents of a file.

        Parameters
        ----------
        filename : str
            The path to the file.

        Returns
        -------
        str
            The checksum of the file contents.
        """
        with open(filename, "r") as f:
            text = f.read()
        checksum = self._text_checksum(text)
        return checksum

    @staticmethod
    def checksum_from_conda_yml_file(self, env_yml, overwrite):
        """
        Generate a checksum for a conda environment YAML file.

        Parameters
        ----------
        env_yml : str
            The path to the conda environment YAML file.
        overwrite : bool
            Whether to overwrite the name and prefix in the YAML file with the checksum.

        Returns
        -------
        str
            The checksum of the YAML file contents.
        """
        with open(env_yml, "r") as f:
            name_idx = None
            pref_idx = None
            R = []
            for i, r in enumerate(f):
                if r[:4] == "name":
                    name_idx = i
                if r[:6] == "prefix":
                    pref_idx = i
                R += [r]
            S = R[(name_idx + 1) : pref_idx]
            text = "".join(S)
        checksum = self._text_checksum(text)
        if overwrite:
            with open(env_yml, "w") as f:
                f.write("name: %s\n" % checksum)
                for s in S:
                    f.write(s)
                prefix = os.path.join(
                    os.sep.join(R[pref_idx].split("prefix: ")[1].split(os.sep)[:-1]),
                    checksum,
                )
                f.write("prefix: %s\n" % prefix)
                f.write("\n")
        return checksum

    def get_conda_and_pip_install_commands_from_dockerfile_if_exclusive(
        self, path, force_exclusive=True
    ):
        """
        Identify install commands from a Dockerfile.

        Parameters
        ----------
        path : str
            The path to the Dockerfile.
        force_exclusive : bool, optional
            Whether to force exclusive conda and pip commands. Default is True.

        Returns
        -------
        list
            A list of install commands.
        """
        dp = SimpleDockerfileParser(path)
        runs = dp.get_runs()
        is_valid = True
        runs_ = []
        for r in runs:
            exec = r.split(" ")[0]
            if exec not in ["conda", "pip", "pip3"]:  # TODO write better
                is_valid = False
            if r.startswith("python -m pip") or r.startswith("python3 -m pip"):
                is_valid = True
            if " -y " not in r and exec == "conda":
                runs_ += [r + " -y"]
            else:
                runs_ += [r]
        if not force_exclusive:
            return runs_
        if is_valid:
            return runs_
        else:
            return None

    def specs_from_dockerfile_as_json(self, dockerfile_dir, dest):
        """
        Write a JSON file with the install requirements inferred from the Dockerfile.

        Parameters
        ----------
        dockerfile_dir : str
            The directory containing the Dockerfile.
        dest : str
            The destination directory for the JSON file.

        Returns
        -------
        str
            The path to the JSON file.
        """
        runs = self.get_conda_and_pip_install_commands_from_dockerfile_if_exclusive(
            dockerfile_dir, force_exclusive=False
        )
        dp = SimpleDockerfileParser(dockerfile_dir)
        bi = dp.baseimage
        sp = bi.split("/")
        if len(sp) == 1:
            org = None
            img = sp[0]
        else:
            org = sp[0]
            img = sp[1]
        if ":" in img:
            sp = img.split(":")
            img = sp[0]
            tag = sp[1]
        else:
            raise Exception("No tag was found in the dockerfile base image.")
        d = defaultdict(list)
        for r in runs:
            result = self._parse_install(r)
            if result is None:
                continue
            k = result["tool"], result["channel"]
            d[k] += result["packages"]
        d = dict((k, sorted(set(v))) for k, v in d.items())
        od = OrderedDict()
        od["base-env"] = self.versions.base_conda_name(org, tag)
        for k in sorted(d.keys()):
            if k[1] is None:
                k_ = k[0]
            else:
                k_ = k[0] + " -c " + k[1]
            od[k_] = d[k]
        json_path = os.path.join(dest, self.SPECS_JSON)
        with open(json_path, "w") as f:
            json.dump(od, f, indent=4)
        return json_path

    def get_base_env(self, path):
        """
        Get the base environment from a JSON file.

        Parameters
        ----------
        path : str
            The path to the JSON file.

        Returns
        -------
        str
            The base environment.
        """
        json_path = os.path.join(path, self.SPECS_JSON)
        with open(json_path, "r") as f:
            data = json.load(f)
        return data["base-env"]

    def checksum_from_dockerfile(self, dockerfile_dir, dest=None):
        """
        Generate a checksum for a Dockerfile.

        Parameters
        ----------
        dockerfile_dir : str
            The directory containing the Dockerfile.
        dest : str, optional
            The destination directory for the checksum file. Default is None.

        Returns
        -------
        str
            The checksum of the Dockerfile.
        """
        if dest is None:
            dest = dockerfile_dir
        filename = os.path.join(dest, self.CHECKSUM_FILE)
        json_path = self.specs_from_dockerfile_as_json(dockerfile_dir, dest=dest)
        checksum = self.checksum_from_file(json_path)
        with open(filename, "w") as f:
            f.write(checksum)
        return checksum

    def specs_from_dockerfile(
        self, dockerfile_dir, dest=None, use_checksum=False, name=None
    ):
        """
        Generate specs from a Dockerfile.

        Parameters
        ----------
        dockerfile_dir : str
            The directory containing the Dockerfile.
        dest : str, optional
            The destination directory for the specs file. Default is None.
        use_checksum : bool, optional
            Whether to use a checksum for the specs. Default is False.
        name : str, optional
            The name to use for the specs. Default is None.

        Returns
        -------
        str
            The name or checksum of the specs.
        """
        if use_checksum:
            return self.checksum_from_dockerfile(dockerfile_dir, dest)  # TODO debug
        else:
            if dest is None:
                dest = dockerfile_dir
            self.specs_from_dockerfile_as_json(
                dockerfile_dir, dest=dest
            )  # TODO remove?
            filename = os.path.join(dest, self.CHECKSUM_FILE)
            with open(filename, "w") as f:
                f.write(name)
            return name

    def activate_base(self):
        """
        Generate a script to activate the base conda environment.

        Returns
        -------
        str
            The script to activate the base conda environment.
        """
        if self.is_base():
            return ""
        snippet = """
        source {0}/etc/profile.d/conda.sh
        conda activate {1}
        """.format(self.conda_prefix(False), BASE)
        return snippet


class SimpleConda(CondaUtils):
    """
    A class to manage conda environments using simple commands.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings in JSON format. Default is None.

    Examples
    --------
    .. code-block:: python

        simple_conda = SimpleConda()
        simple_conda.create("myenv", "3.8")

    """

    def __init__(self, config_json=None):
        CondaUtils.__init__(self, config_json=config_json)

    def _env_list(self):
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "env_list.tsv")
        tmp_script = os.path.join(tmp_folder, "script.sh")
        bash_script = """
        source {0}/etc/profile.d/conda.sh
        conda env list > {1}
        """.format(self.conda_prefix(self.is_base()), tmp_file)
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        run_command("bash {0}".format(tmp_script))
        with open(tmp_file, "r") as f:
            envs = []
            for l in f:
                envs += [l.rstrip()]
        return envs

    def create(self, environment, python_version):
        """
        Create a new conda environment.

        Parameters
        ----------
        environment : str
            The name of the new environment.
        python_version : str
            The Python version to use in the new environment.

        Returns
        -------
        None
        """
        cmd = "conda create -n {0} python={1} -y".format(environment, python_version)
        run_command(cmd)

    def active_env(self):
        """
        Get the currently active conda environment.

        Returns
        -------
        str
            The name of the currently active conda environment.
        """
        envs = self._env_list()
        for l in envs:
            if "*" in l:
                return l.split()[0]
        return None

    def exists(self, environment):
        """
        Check if a conda environment exists.

        Parameters
        ----------
        environment : str
            The name of the environment to check.

        Returns
        -------
        bool
            True if the environment exists, False otherwise.
        """
        envs = self._env_list()
        n = len(environment)
        for l in envs:
            if l[:n] == environment and l[n] == " ":
                return True
        return False

    def startswith(self, environment):
        """
        Get a list of conda environments that start with the given prefix.

        Parameters
        ----------
        environment : str
            The prefix to check.

        Returns
        -------
        list
            A list of conda environments that start with the given prefix.
        """
        envs = self._env_list()
        envs_list = []
        for l in envs:
            if l.startswith(environment):
                envs_list += [l.split(" ")[0]]
        return envs_list

    def get_python_path_env(self, environment):
        """
        Get the Python path for a conda environment.

        Parameters
        ----------
        environment : str
            The name of the environment.

        Returns
        -------
        str
            The Python path for the environment.
        """
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "tmp.txt")
        self.run_commandlines(environment, "which python > {0}".format(tmp_file))
        with open(tmp_file, "r") as f:
            python_path = f.read().rstrip()
        return python_path

    def delete_one(self, environment):
        """
        Delete a conda environment.

        Parameters
        ----------
        environment : str
            The name of the environment to delete.
        """
        if not self.exists(environment):
            return
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_script = os.path.join(tmp_folder, "script.sh")
        bash_script = self.activate_base()
        bash_script += """
        source {0}/etc/profile.d/conda.sh
        conda env remove --name {1} -y
        """.format(self.conda_prefix(True), environment)
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        run_command("bash {0}".format(tmp_script))

    def _proper_delete(self, environment):
        for env in self.startswith(environment):
            self.delete_one(env)

    def _manual_delete(self, environment):
        envs_path = os.path.join(self.conda_prefix(True), "envs")
        for env in os.listdir(envs_path):
            if env.startswith(environment):
                shutil.rmtree(os.path.join(envs_path, env))

    def delete(self, environment):
        """
        Delete a conda environment.

        Parameters
        ----------
        environment : str
            The name of the environment to delete.
        """
        self._proper_delete(environment)
        self._manual_delete(environment)

    def export_env_yml(self, environment, dest):
        """
        Export a conda environment as an environment YAML file.

        Parameters
        ----------
        environment : str
            The name of the environment to export.
        dest : str
            The destination directory for the YAML file.
        """
        if not self.exists(environment):
            return
        if self.is_base():
            return
        yml_file = os.path.join(dest, CONDA_ENV_YML_FILE)
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_script = os.path.join(tmp_folder, "script.sh")
        bash_script = self.activate_base()
        bash_script += """
        source {0}/etc/profile.d/conda.sh
        conda activate {1}
        conda env export --no-builds > {2}
        conda deactivate
        """.format(self.conda_prefix(True), environment, yml_file)
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        run_command("bash {0}".format(tmp_script))

    def clone(self, src_env, dst_env):
        """
        Make an exact copy of a conda environment.

        Parameters
        ----------
        src_env : str
            The name of the source environment.
        dst_env : str
            The name of the destination environment.

        Raises
        ------
        Exception
            If the source environment does not exist or the destination environment already exists.
        """
        if not self.exists(src_env):
            raise Exception("{0} source environment does not exist".format(src_env))
        if self.exists(dst_env):
            raise Exception("{0} destination environment exists".format(dst_env))
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_script = os.path.join(tmp_folder, "script.sh")
        bash_script = self.activate_base()
        bash_script += """
        source {0}/etc/profile.d/conda.sh
        conda create --clone {1} --name {2} -y
        """.format(self.conda_prefix(True), src_env, dst_env)
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        run_command("bash {0}".format(tmp_script))

    def _catch_critical_errors_in_conda(self, log):
        log = log.split(os.linesep)
        critical_errors = []
        for l in log:
            if l.startswith("ERROR"):
                if "dependency resolver" not in l:
                    critical_errors += [l]
        return critical_errors

    def create_executable_bash_script(self, environment, commandlines, file_name):
        """
        Create an executable bash script to run commands in a conda environment.

        Parameters
        ----------
        environment : str
            The name of the environment.
        commandlines : str or list
            The commands to run.
        file_name : str
            The name of the bash script file.

        Returns
        -------
        str
            The path to the bash script file.
        """
        if type(commandlines) is list:
            commandlines = "\n".join(commandlines)
        bash_script = self.activate_base()
        bash_script += """
        source {0}/etc/profile.d/conda.sh
        conda activate {1}
        {2}
        """.format(self.conda_prefix(True), environment, commandlines)
        with open(file_name, "w") as f:
            f.write(bash_script)
        return file_name

    @throw_ersilia_exception()
    def run_commandlines(self, environment, commandlines):
        """
        Run commands in a given conda environment.

        Parameters
        ----------
        environment : str
            The name of the environment.
        commandlines : str or list
            The commands to run.

        Raises
        ------
        ModelPackageInstallError
            If critical errors are found in the conda logs.
        """
        if type(commandlines) is list:
            commandlines = " && ".join(commandlines)
        logger.debug("Run commandlines on {0}".format(environment))
        logger.debug(commandlines)
        if not self.exists(environment):
            raise Exception("{0} environment does not exist".format(environment))
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_script = os.path.join(tmp_folder, "script.sh")
        logger.debug("Activating base environment")
        logger.debug("Current working directory: {0}".format(os.getcwd()))
        self.create_executable_bash_script(environment, commandlines, tmp_script)
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_log = os.path.join(tmp_folder, "command_outputs.log")
        cmd = "bash {0} 2>&1 | tee -a {1}".format(tmp_script, tmp_log)
        logger.debug("Running {0}".format(cmd))
        run_command(cmd)
        with open(tmp_log, "r") as f:
            log_file = f.read()
        logger.debug(log_file)
        if "ERROR" in log_file:
            logger.debug("Error occurred while running: {0}".format(cmd))
            critical_errors = self._catch_critical_errors_in_conda(log_file)
            if len(critical_errors) > 0:
                raise ModelPackageInstallError(cmd)
        logger.debug("Activation done")


class StandaloneConda(object):
    """
    A class to manage standalone conda environments.
    """

    def __init__(self):
        pass

    def exists(self, environment):
        """
        Check if a standalone conda environment exists.

        Parameters
        ----------
        environment : str
            The name of the environment to check.

        Returns
        -------
        bool
            True if the environment exists, False otherwise.
        """
        return os.path.exists(os.path.join("/", environment))

    def run_commandlines(self, environment, commandlines):
        """
        Run commands in a given standalone conda environment.

        Parameters
        ----------
        environment : str
            The name of the environment.
        commandlines : str or list
            The commands to run.
        """
        logger.debug("Run commandlines on {0}".format(environment))
        logger.debug(commandlines)

        if not self.exists(environment):
            raise Exception("{0} environment does not exist".format(environment))

        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_script = os.path.join(tmp_folder, "script.sh")
        logger.info(f"Activating environment: {environment}")
        logger.debug("Current working directory: {0}".format(os.getcwd()))
        bash_script = """
        source /{0}/bin/activate
        {1}
        """.format(environment, commandlines)
        with open(tmp_script, "w") as f:
            f.write(bash_script)

        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_log = os.path.join(tmp_folder, "command_outputs.log")
        cmd = "bash {0} 2>&1 | tee -a {1}".format(tmp_script, tmp_log)
        logger.debug("Running {0}".format(cmd))
        run_command(cmd)
        with open(tmp_log, "r") as f:
            log_file = f.read()
        logger.debug(log_file)
        if "ERROR" in log_file:
            logger.debug("Error occurred while running: {0}".format(cmd))
            # TODO do we catch critical errors here?
        logger.debug("Activation done")
