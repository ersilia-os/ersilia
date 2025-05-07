import os
import subprocess
import sys
import warnings
import zipfile
from typing import Any
import warnings

warnings.filterwarnings("ignore", message="Using slow pure-python SequenceMatcher")

from ....default import (
    EOS_TMP,
    GITHUB_ORG,
    S3_BUCKET_URL_ZIP,
)
from ....hub.content.card import ModelCard
from ....utils.download import GitHubDownloader, S3Downloader
from ....utils.conda import SimpleConda
from ....utils.logging import make_temp_dir
from ....utils.terminal import yes_no_input
from ....cli import echo

warnings.filterwarnings("ignore", message="Using slow pure-python SequenceMatcher.*")


class SetupService:
    """
    Service for setting up the environment and fetching the model repository.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    dir : str
        Directory where the model repository will be cloned.
    from_github : bool
        Flag indicating whether to fetch the repository from GitHub.
    from_s3 : bool
        Flag indicating whether to fetch the repository from S3.
    logger : Any
        Logger for logging messages.
    """

    BASE_URL = "https://github.com/ersilia-os/"

    def __init__(
        self,
        model_id: str,
        dir: str,
        from_github: bool,
        from_s3: bool,
        logger: Any,
    ):
        self.model_id = model_id
        self.dir = dir
        self.logger = logger
        self.from_github = from_github
        self.from_s3 = from_s3

        self.mc = ModelCard()
        self.metadata = self.mc.get(model_id)
        self.s3 = self.metadata.get("card", {}).get("S3") or self.metadata.get(
            "metadata", {}
        ).get("S3")
        self.repo_url = f"{self.BASE_URL}{self.model_id}"
        self.overwrite = self._handle_overwrite()
        self.github_down = GitHubDownloader(overwrite=self.overwrite)
        self.s3_down = S3Downloader()
        self.conda = SimpleConda()

    def _handle_overwrite(self) -> bool:
        if os.path.exists(self.dir):
            self.logger.info(f"Directory {self.dir} already exists.")
            return yes_no_input(
                f"Directory {self.dir} already exists. Do you want to overwrite it? [Y/n]",
                default_answer="n",
            )
        return False

    def _download_s3(self):
        if not self.overwrite and os.path.exists(self.dir):
            self.logger.info("Skipping S3 download as user chose not to overwrite.")
            return

        tmp_file = os.path.join(make_temp_dir("ersilia-"), f"{self.model_id}.zip")

        self.logger.info(f"Downloading model from S3 to temporary file: {tmp_file}")
        self.s3_down.download_from_s3(
            bucket_url=S3_BUCKET_URL_ZIP,
            file_name=f"{self.model_id}.zip",
            destination=tmp_file,
        )

        self.logger.info(f"Extracting model to: {self.dir}")
        with zipfile.ZipFile(tmp_file, "r") as zip_ref:
            zip_ref.extractall(EOS_TMP)

    def _download_github(self):
        try:
            if not os.path.exists(EOS_TMP):
                self.logger.info(f"Path does not exist. Creating: {EOS_TMP}")
                os.makedirs(EOS_TMP, exist_ok=True)
        except OSError as e:
            self.logger.error(f"Failed to create directory {EOS_TMP}: {e}")

        self.logger.info(f"Cloning repository from GitHub to: {EOS_TMP}")
        self.github_down.clone(
            org=GITHUB_ORG,
            repo=self.model_id,
            destination=self.dir,
        )

    def get_model(self):
        if self.from_s3:
            self._download_s3()

        if self.from_github:
            self._download_github()

    @staticmethod
    def run_command(
        command: str, logger, capture_output: bool = False, shell: bool = True
    ) -> str:
        """
        Run a shell command.

        Parameters
        ----------
        command : str
            The command to run.
        logger : logging.Logger
            Logger for logging messages.
        capture_output : bool, optional
            Flag indicating whether to capture the command output.
        shell : bool, optional
            Flag indicating whether to run the command in the shell.

        Returns
        -------
        str
            The output of the command.

        Raises
        ------
        subprocess.CalledProcessError
            If the command returns a non-zero exit code.
        """
        try:
            if capture_output:
                result = subprocess.run(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                    shell=shell,
                )
                return result.stdout
            else:
                process = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    shell=shell,
                )

                stdout_lines, stderr_lines = [], []

                for line in iter(process.stdout.readline, ""):
                    if line.strip():
                        stdout_lines.append(line.strip())
                        logger.info(line.strip())

                for line in iter(process.stderr.readline, ""):
                    if line.strip():
                        stderr_lines.append(line.strip())
                        logger.info(line.strip())

                process.wait()
                return process.stdout

        except subprocess.CalledProcessError as e:
            logger.error(f"Error executing command: {e}")
            echo(f"\n Error happened when : {e}", fg="red")
            if e.output:
                logger.debug(f"Output: {e.output.strip()}")
                echo(f"Output: {e.output.strip()}", fg="red")
            if e.stderr:
                logger.error(f"Error: {e.stderr.strip()}")
                echo(f"Error: {e.stderr.strip()}", fg="red")
            sys.exit(1)
        except Exception as e:
            logger.debug(f"Unexpected error: {e}")
            echo(f"\n Unexpected error: {e}", fg="red")
            sys.exit(1)

    @staticmethod
    def get_conda_env_location(model_id: str, logger) -> str:
        """
        Get the location of the Conda environment for the model.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        logger : logging.Logger
            Logger for logging messages.

        Returns
        -------
        str
            The location of the Conda environment.

        Raises
        ------
        subprocess.CalledProcessError
            If the command to list Conda environments returns a non-zero exit code.
        """
        try:
            result = SetupService.run_command(
                "conda env list", logger=logger, capture_output=True
            )
            for line in result.splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if parts[0] == model_id:
                    return parts[-1]
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running conda command: {e.stderr}")
            echo(f"Error running conda command: {e.stderr}", fg="red")
        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            echo(f"Unexpected error: {e}", fg="red")

        return None
