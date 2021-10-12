import os

try:
    from isaura.default import REPOSITORY_PATH as ISAURA_REPOSITORY_PATH
except:
    ISAURA_REPOSITORY_PATH = None

from .. import ErsiliaBase


class LakeBase(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        if ISAURA_REPOSITORY_PATH is not None:
            self.lake_dir = os.path.abspath(ISAURA_REPOSITORY_PATH)
        else:
            self.lake_dir = None
