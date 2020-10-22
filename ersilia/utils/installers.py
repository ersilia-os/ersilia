import shutil
import subprocess
import os


class Installer(object):

    def __init__(self):
        pass

    @staticmethod
    def _is_tool(name):
        return shutil.which(name) is not None

    def conda(self):
        if self._is_tool("conda"):
            return
        subprocess.Popen("pip install -y conda", shell=True).wait()

    def git(self):
        if self._is_tool("git"):
            return
        self.conda()
        subprocess.Popen("conda install -y -q git", shell=True).wait()

    def rdkit(self):
        try:
            import rdkit
            exists = True
        except ModuleNotFoundError:
            exists = False
        if exists:
            return
        subprocess.Popen("conda install -c conda-forge -y -q rdkit", shell=True).wait()

    def config(self):
        from ..default import EOS
        if os.path.exists(os.path.join(EOS, "config.json")):
            return
        os.makedirs(EOS, exist_ok=True)
        from .download import GitHubDownloader
        gd = GitHubDownloader(overwrite=True)
        gd.download_single("ersilia-os", "ersilia", "config.json", os.path.join(EOS, "config.json"))


def check_dependencies():
    ins = Installer()
    ins.conda()
    ins.git()
    ins.rdkit()
    ins.config()
