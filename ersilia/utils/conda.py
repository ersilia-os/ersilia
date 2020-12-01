import tempfile
import os
import json
from collections import defaultdict, OrderedDict
from .terminal import run_command
from .docker import SimpleDockerfileParser


BASE = "base"


class BaseConda(object):

    def __init__(self):
        pass

    @staticmethod
    def default_env():
        return os.environ["CONDA_DEFAULT_ENV"]

    def is_base(self):
        default_env = self.default_env()
        if default_env == BASE:
            return True
        else:
            return False

    @staticmethod
    def conda_prefix(is_base):
        if is_base:
            return "CONDA_PREFIX"
        else:
            return "CONDA_PREFIX_1"


class CondaBashSnippets(BaseConda):

    def __init__(self):
        BaseConda.__init__(self)

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
        result = {
            "tool": tool,
            "channel": channel,
            "packages": packages
        }
        return result

    def install_from_dockerfile(self, path):
        """Identifies install commands from Dockerfile
        Right now this command is conservative and returns None if at least
        one of the commands is not conda ... or pip ... or pip3 ...
        """
        dp = SimpleDockerfileParser(path)
        runs = dp.get_runs()
        is_valid = True
        for r in runs:
            exec = r.split(" ")[0]
            if exec not in ["conda", "pip", "pip3"]:
                is_valid = False
        if is_valid:
            return runs
        else:
            return None

    def specs_from_dockerfile_as_json(self, dockerfile_path, dest):
        """Writes a json file with the install requirements inferred from the Dockerfile.
        """
        runs = self.install_from_dockerfile(dockerfile_path)
        if not runs:
            return None
        dp = SimpleDockerfileParser(dockerfile_path)
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
            tag = None
        d = defaultdict(list)
        for r in runs:
            result = self._parse_install(r)
            if result is None:
                continue
            k = result["tool"], result["channel"]
            d[k] += result["packages"]
        d = dict((k, sorted(set(v))) for k,v in d.items())
        od = OrderedDict()
        od["basesenv"] = img
        for k in sorted(d.keys()):
            if k[1] is None:
                k_ = k[0]
            else:
                k_ = k[0] + " -c " + k[1]
            od[k_] = d[k]
        json_path = os.path.join(dest, "specs.json")
        with open(json_path, "w") as f:
            json.dump(od, f, indent=4)

    def activate_base(self):
        if self.is_base():
            return ""
        snippet = """
        source ${0}/etc/profile.d/conda.sh
        conda activate {1}
        """.format(
            self.conda_prefix(False),
            BASE
        )
        return snippet


class SimpleConda(CondaBashSnippets):

    def __init__(self):
        CondaBashSnippets.__init__(self)

    def _env_list(self):
        tmp_folder = tempfile.mkdtemp()
        tmp_file = os.path.join(tmp_folder, "env_list.tsv")
        tmp_script = os.path.join(tmp_folder, "script.sh")
        bash_script = """
        source ${0}/etc/profile.d/conda.sh
        conda env list > {1}
        """.format(
            self.conda_prefix(self.is_base()),
            tmp_file
        )
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        run_command("bash {0}".format(tmp_script), quiet=True)
        with open(tmp_file, "r") as f:
            envs = []
            for l in f:
                envs += [l.rstrip()]
        return envs

    def active_env(self):
        envs = self._env_list()
        for l in envs:
            if "*" in l:
                return l.split()[0]
        return None

    def exists(self, environment):
        envs = self._env_list()
        n = len(environment)
        for l in envs:
            if l[:n] == environment:
                return True
        return False

    def export_env_yml(self, environment, dest):
        """
        Export conda environment as an environment.yml file.
        """
        if not self.exists(environment):
            return
        if self.is_base():
            return
        yml_file = os.path.join(dest, "environment.yml")
        tmp_folder = tempfile.mkdtemp()
        tmp_script = os.path.join(tmp_folder, "script.sh")
        bash_script = self.activate_base()
        bash_script += """
        source ${0}/etc/profile.d/conda.sh
        conda activate {1}
        conda env export > {2}
        conda deactivate
        """.format(
            self.conda_prefix(True),
            environment,
            yml_file
        )
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        run_command("bash {0}".format(tmp_script), quiet=True)
