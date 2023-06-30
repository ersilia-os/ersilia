import os

from .. import ErsiliaBase, EOS
from ..core.session import ERSILIA_RUNS_FOLDER


class S3Logger(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.runs_directory = os.path.join(EOS, ERSILIA_RUNS_FOLDER)

    def upload(self):
        pass
