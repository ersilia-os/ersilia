import os
import boto3

from .. import ErsiliaBase, EOS
from ..core.session import ERSILIA_RUNS_FOLDER


class S3Logger(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.runs_directory = os.path.join(EOS, ERSILIA_RUNS_FOLDER)
        self.aws_access_key_id = None
        self.aws_secret_access_key = None

    def set_credentials(self, aws_access_key_id, aws_secret_access_key):
        self.aws_access_key_id = aws_access_key_id
        self.aws_secret_access_key = aws_secret_access_key

    def _upload_log(self):
        pass

    def _upload_lake(self):
        pass

    def _upload_meta(self):
        pass

    def upload(self):
        if self.aws_access_key_id is None or self.aws_secret_access_key is None:
            self.logger.warning(
                "It was not possible to upload to S3. AWS access key or secret access key was not provided. Please use the set_credentials method."
            )
        self._upload_log()
        self._upload_lake()
        self._upload_meta()
