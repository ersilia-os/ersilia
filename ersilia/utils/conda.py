import tempfile
import os
from .terminal import run_command


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




class SimpleConda(BaseConda):

    def __init__(self):
        BaseConda.__init__(self)
        self.snips = CondaBashSnippets()

    def exists(self, environment):
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
            n = len(environment)
            for l in f:
                if l[:n] == environment:
                    return True
        return False
