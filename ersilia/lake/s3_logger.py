import os

from .. import EOS, ErsiliaBase
from ..default import ERSILIA_RUNS_FOLDER


class S3Logger(ErsiliaBase):
    """
    Logger class for uploading logs to AWS S3.

    Parameters
    ----------
    model_id : str
        Identifier for the model.
    config_json : dict, optional
        Configuration settings in JSON format.

    Attributes
    ----------
    model_id : str
        Identifier for the model.
    runs_directory : str
        Directory path for storing run logs.
    aws_access_key_id : str or None
        AWS access key ID.
    aws_secret_access_key : str or None
        AWS secret access key.
    """

    def __init__(self, model_id: str, config_json: dict = None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.runs_directory = os.path.join(EOS, ERSILIA_RUNS_FOLDER)
        self.aws_access_key_id = None
        self.aws_secret_access_key = None

    def set_credentials(self, aws_access_key_id: str, aws_secret_access_key: str):
        """
        Set AWS credentials for S3 access.

        Parameters
        ----------
        aws_access_key_id : str
            AWS access key ID.
        aws_secret_access_key : str
            AWS secret access key.
        """
        self.aws_access_key_id = aws_access_key_id
        self.aws_secret_access_key = aws_secret_access_key

    def _upload_log(self):
        pass

    def _upload_lake(self):
        pass

    def _upload_meta(self):
        pass

    def upload(self):
        """
        Upload logs, lake files, and metadata to S3.

        Warns
        -----
        Warning
            If AWS credentials are not set, a warning is logged.
        """
        if self.aws_access_key_id is None or self.aws_secret_access_key is None:
            self.logger.warning(
                "It was not possible to upload to S3. AWS access key or secret access key was not provided. Please use the set_credentials method."
            )
        self._upload_log()
        self._upload_lake()
        self._upload_meta()
