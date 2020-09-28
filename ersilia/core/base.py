from ..utils.config import Config, Credentials
import os


class ErsiliaBase(object):

    def __init__(self, config_json=None, credentials_json=None):
        self.cfg = Config(json_file=config_json)
        self.cred = Credentials(json_file=credentials_json)
        self._tmp_dir = self._abs_path(self.cfg.LOCAL.TMP)
        if not os.path.exists(self._tmp_dir):
            os.makedirs(self._tmp_dir, exist_ok=True)
        self._dest_dir = self._abs_path(self.cfg.LOCAL.DEST)
        if not os.path.exists(self._dest_dir):
            os.makedirs(self._dest_dir, exist_ok=True)

    @staticmethod
    def _abs_path(path):
        if path[0] == "~":
            home = os.path.expanduser("~")
            return home+path[1:]
        return os.path.abspath(path)
