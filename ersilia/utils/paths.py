import re
import os
import collections
from pathlib import Path

MODELS_DEVEL_DIRNAME = "models"


class Paths(object):
    def __init__(self):
        self.essentials = ["setup.py", "README.md", "CODE_OF_CONDUCT.md"]

    @staticmethod
    def _eos_regex():
        return re.compile(r"eos[0-9][a-z0-9]{3}")

    @staticmethod
    def home():
        """Get home directory"""
        return os.path.abspath(str(Path.home()))

    def model_id_from_path(self, path):
        """Guess model identifier based on the path"""
        regex = self._eos_regex()
        path = os.path.abspath(path)
        model_ids = sorted(set(regex.findall(path)))
        if len(model_ids) == 1:
            return model_ids[0]
        else:
            return None

    def org_development_path(self):
        """Guess generic development path"""
        path = self.ersilia_development_path()
        if path is None:
            return
        else:
            return os.path.split(path)[0]

    def ersilia_development_path(self):
        """Try to guess the package development path in the local computer"""
        path = os.path.dirname(__file__)
        for _ in range(2):
            path = os.path.split(path)[0]
        for essential in self.essentials:
            if not os.path.exists(os.path.join(path, essential)):
                return None
        return path

    def models_development_path(self):
        """Try to guess the models path in the local computer.
        The directory with more 'eos0xxx' subdirectories will be returned."""
        org_dev_path = self.org_development_path()
        if org_dev_path is None:
            return
        regex = self._eos_regex()
        cands = collections.defaultdict(int)
        for dirpath, dirnames, filenames in os.walk(org_dev_path):
            ap = os.path.abspath(dirpath)
            bn = os.path.basename(ap)
            if bn == MODELS_DEVEL_DIRNAME:
                for dn in dirnames:
                    if regex.match(dn):
                        cands[ap] += 1
        path = sorted(cands.items(), key=lambda item: -item[1])[0][0]
        return path

    @staticmethod
    def exists(path):
        if path is None:
            return False
        if os.path.exists(path):
            return True
        else:
            return False
