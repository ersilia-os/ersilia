import os
import subprocess
from pathlib import Path
from ..utils.config import Config, Credentials
from ..default import EOS
from .. import logger

home = str(Path.home())


class ErsiliaBase(object):
    """Base class of Ersilia.

    This class is used as a configuration for many of the classes of the package.
    """

    def __init__(self, config_json=None, credentials_json=None):
        """Initialize Ersilia base class.

        Attributes:
            config_json (str, optional): Path to a configuration file.
            credentials_json (str, optional): Path to a credentials file.
        """
        self.eos_dir = EOS
        self.config_json = config_json
        self.credentials_json = credentials_json
        self.cfg = Config(json_file=config_json)
        self.cred = Credentials(json_file=credentials_json)
        self._lake_dir = self._abs_path(os.path.join(EOS, "isaura", "lake"))
        self._tmp_dir = self._abs_path(os.path.join(EOS, self.cfg.LOCAL.TMP))
        if not os.path.exists(self._tmp_dir):
            os.makedirs(self._tmp_dir, exist_ok=True)
        self._dest_dir = self._abs_path(os.path.join(EOS, self.cfg.LOCAL.DEST))
        if not os.path.exists(self._dest_dir):
            os.makedirs(self._dest_dir, exist_ok=True)
        self._bentoml_dir = os.path.join(
            self._abs_path(os.path.join(Path.home(), "bentoml")), "repository"
        )
        self._bundles_dir = os.path.join(self.eos_dir, "repository")
        if not os.path.exists(self._bundles_dir):
            os.makedirs(self._bundles_dir, exist_ok=True)
        self.logger = logger

    @staticmethod
    def _abs_path(path):
        return os.path.abspath(path)

    def _model_path(self, model_id):
        folder = os.path.join(self._dest_dir, model_id)
        return folder

    def _get_latest_bentoml_tag(self, model_id):
        path = os.path.join(self._bentoml_dir, model_id)
        if not os.path.exists(path):
            return None
        items = sorted(os.listdir(path))
        if not items:
            return None
        else:
            return items[-1]

    def _get_latest_bundle_tag(self, model_id):
        path = os.path.join(self._bundles_dir, model_id)
        if not os.path.exists(path):
            return None
        items = sorted(os.listdir(path))
        if not items:
            return None
        else:
            return items[-1]

    def _get_bentoml_location(self, model_id):
        tag = self._get_latest_bentoml_tag(model_id)
        path = os.path.join(self._bentoml_dir, model_id)
        if not os.path.exists(path):
            return None
        if tag is not None:
            return os.path.join(path, tag)
        else:
            return path

    def _get_bundle_location(self, model_id):
        tag = self._get_latest_bundle_tag(model_id)
        path = os.path.join(self._bundles_dir, model_id)
        if not os.path.exists(path):
            return None
        if tag is not None:
            return os.path.join(path, tag)
        else:
            return path

    @staticmethod
    def _get_bento_location(model_id):
        cmd = ["bentoml", "get", "%s:latest" % model_id, "--print-location", "--quiet"]
        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        result = result.stdout.decode("utf-8").rstrip()
        return result

    def _is_ready(self, model_id):
        """Check whether a model exists in the local computer"""
        try:
            self._get_latest_bundle_tag(model_id)
        except:
            return False
        path = os.path.join(self._abs_path(self._dest_dir), model_id)
        if not os.path.exists(path):
            return False
        return True

    def _has_credentials(self):
        if self.cred is None:
            return False
        return True
