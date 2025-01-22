import os
import shutil
import zipfile

import boto3

from .. import ErsiliaBase
from ..default import ERSILIA_MODELS_S3_BUCKET, ERSILIA_MODELS_ZIP_S3_BUCKET
from ..utils.logging import make_temp_dir
from ..utils.terminal import run_command

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

        uploader = S3BucketRepoUploader(
            model_id="model_id",
            config_json="path/to/config.json",
        )
        uploader.set_credentials(
            aws_access_key_id="access_key",
            aws_secret_access_key="secret_key",
        )
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

    def _clone(self):
        self.logger.debug("Cloning model {0} from ersilia-os".format(self.model_id))
        run_command(
            "cd {0}; git clone https://github.com/ersilia-os/{1}; cd {2}".format(
                self.tmp_folder, self.model_id, self.cwd
            )
        )

    def _ungit(self, repo_path):
        self.logger.debug("Removing git files")
        dotgit_folder = os.path.join(repo_path, ".git")
        gitignore_file = os.path.join(repo_path, ".gitignore")
        github_file = os.path.join(repo_path, ".github")
        if os.path.exists(dotgit_folder):
            shutil.rmtree(dotgit_folder)
        if os.path.exists(gitignore_file):
            os.remove(gitignore_file)
        if os.path.exists(github_file):
            shutil.rmtree(github_file)

    def _upload_files(self, repo_path):
        self.logger.debug("Uploading repo files")
        repo_path = os.path.abspath(repo_path)
        basename = os.path.basename(repo_path)
        self.logger.debug("Taking basename {0}".format(basename))
        session = boto3.Session(
            aws_access_key_id=self.aws_access_key_id,
            aws_secret_access_key=self.aws_secret_access_key,
            region_name=AWS_ACCOUNT_REGION,
        )
        self.logger.debug(session)
        s3 = session.resource("s3")
        bucket = s3.Bucket(ERSILIA_MODELS_S3_BUCKET)
        self._delete_model_from_s3(bucket)
        model_id = self.model_id
        for subdir, _, files in os.walk(repo_path):
            for file in files:
                if file in self.ignore:
                    self.logger.debug("Ignoring {0}".format(file))
                    continue
                full_path = os.path.join(subdir, file)
                with open(full_path, "rb") as data:
                    s = full_path.split(basename)[-1]
                    if not s.startswith("/"):
                        s = "/" + s
                    self.logger.debug(s)
                    key = model_id + s
                    self.logger.debug(key)
                    bucket.put_object(Key=key, Body=data, ACL="public-read")

    def _delete_model_from_s3(self, bucket):
        bucket.objects.filter(Prefix="{0}/".format(self.model_id)).delete()

    def _zipdir(self, repo_path, ziph):
        for root, dirs, files in os.walk(repo_path):
            for file in files:
                if file in self.ignore:
                    continue
                ziph.write(
                    os.path.join(root, file),
                    os.path.relpath(
                        os.path.join(root, file), os.path.join(repo_path, "..")
                    ),
                )

    def _zip_model(self, repo_path):
        repo_path = os.path.abspath(repo_path)
        self.zip_model_file = os.path.join(self.tmp_zip_folder, self.model_id + ".zip")
        zipf = zipfile.ZipFile(self.zip_model_file, "w", zipfile.ZIP_DEFLATED)
        self._zipdir(repo_path, zipf)
        zipf.close()

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
