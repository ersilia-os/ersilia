import boto3
import os
import shutil
import tempfile
import zipfile

from ..utils.terminal import run_command
from ..utils.logging import make_temp_dir
from .. import ErsiliaBase
from ..default import ERSILIA_MODELS_S3_BUCKET, ERSILIA_MODELS_ZIP_S3_BUCKET

AWS_ACCOUNT_REGION = "eu-central-1"


class S3BucketRepoUploader(ErsiliaBase):
    """
    Class for uploading model repositories to an S3 bucket.

    Parameters
    ----------
    model_id : str
        The ID of the model to be uploaded.
    config_json : str, optional
        Path to the configuration JSON file.

    Examples
    --------
    .. code-block:: python

        uploader = S3BucketRepoUploader(model_id="model_id", config_json="path/to/config.json")
        uploader.set_credentials(aws_access_key_id="access_key", aws_secret_access_key="secret_key")
        uploader.upload()
    """
    def __init__(self, model_id: str, config_json=None):
        self.model_id = model_id
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.cwd = os.getcwd()
        self.tmp_folder = make_temp_dir(prefix="ersilia-")
        self.tmp_zip_folder = make_temp_dir(prefix="ersilia-")
        self.aws_access_key_id = None
        self.aws_secret_access_key = None
        self.ignore = ["upload_model_to_s3.py"]

    def set_credentials(self, aws_access_key_id: str, aws_secret_access_key: str):
        """
        Set AWS credentials.

        Parameters
        ----------
        aws_access_key_id : str
            AWS access key ID.
        aws_secret_access_key : str
            AWS secret access key.
        """
        self.aws_access_key_id = aws_access_key_id
        self.aws_secret_access_key = aws_secret_access_key

    def upload(self, repo_path=None):
        """
        Upload the model repository to the S3 bucket.

        Parameters
        ----------
        repo_path : str, optional
            Path to the local repository. If not provided, the repository will be cloned from GitHub.
        """
        if repo_path is not None:
            self.logger.debug("Repo path is {0}".format(os.path.abspath(repo_path)))
            self._ungit(repo_path=repo_path)
        else:
            repo_path = os.path.join(self.tmp_folder, self.model_id)
            self._clone()
            self._ungit(repo_path=repo_path)
        self.logger.debug(
            "Uploading model folder to S3 bucket {0}".format(ERSILIA_MODELS_S3_BUCKET)
        )
        self._upload_files(repo_path)

    def upload_zip(self, repo_path=None):
        """
        Upload the zipped model repository to the S3 bucket.

        Parameters
        ----------
        repo_path : str, optional
            Path to the local repository. If not provided, the repository will be cloned from GitHub.
        """
        if repo_path is not None:
            self.logger.debug("Repo path is {0}".format(os.path.abspath(repo_path)))
            self._ungit(repo_path=repo_path)
        else:
            repo_path = os.path.join(self.tmp_folder, self.model_id)
            self._clone()
            self._ungit(repo_path=repo_path)
        self.logger.debug(
            "Uploading zipped model folder to S3 bucket {0}".format(
                ERSILIA_MODELS_S3_BUCKET
            )
        )
        self._zip_model(repo_path)
        session = boto3.Session(
            aws_access_key_id=self.aws_access_key_id,
            aws_secret_access_key=self.aws_secret_access_key,
            region_name=AWS_ACCOUNT_REGION,
        )
        self.logger.debug(session)
        s3 = session.client("s3")
        key = self.model_id
        self.logger.debug(key)
        s3.upload_file(
            self.zip_model_file,
            ERSILIA_MODELS_ZIP_S3_BUCKET,
            self.model_id + ".zip",
            ExtraArgs={"ACL": "public-read"},
        )
