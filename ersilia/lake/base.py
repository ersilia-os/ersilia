import os

from isaura.default import REPOSITORY_PATH as ISAURA_REPOSITORY_PATH
from .. import ErsiliaBase


class LakeBase(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.lake_dir = os.path.abspath(ISAURA_REPOSITORY_PATH)
