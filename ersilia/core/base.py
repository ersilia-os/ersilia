from ..utils.config import Config, Credentials
import os


class ErsiliaBase(object):

    def __init__(self, config_json=None, credentials_json=None):
        self.config_json = config_json
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

    def _get_latest_bentoml_tag(self, model_id):
        path = os.path.join(self._abs_path(self.cfg.LOCAL.BENTOML), "repository", model_id)
        return sorted(os.listdir(path))[-1]

    def _is_ready(self, model_id):
        """Check whether a model exists in the local computer"""
        try:
            self._get_latest_bentoml_tag(model_id)
        except:
            return False
        path = os.path.join(self._abs_path(self.cfg.LOCAL.DEST), model_id)
        if not os.path.exists(path):
            return False
        return True

    def _has_credentials(self):
        if self.cred is None:
            return False
        return True
