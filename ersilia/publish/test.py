import csv
import docker
import json
import os
import subprocess
import zipfile
import sys
import tempfile
import math
from pathlib import Path
import yaml
import sys
import types
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from typing import Callable
from typing import List, Any

# ruff: noqa
MISSING_PACKAGES = False
try:
    from fuzzywuzzy import fuzz
    from rich.console import Console
    from rich.table import Table
    from rich.text import Text
    from scipy.stats import spearmanr
except ImportError:
    MISSING_PACKAGES = True
# ruff: enable

from .. import ErsiliaBase, throw_ersilia_exception
from ..default import (
    DOCKERFILE_FILE,
    INSTALL_YAML_FILE,
    METADATA_JSON_FILE,
    METADATA_YAML_FILE,
    PACK_METHOD_BENTOML,
    PACK_METHOD_FASTAPI,
    PREDEFINED_EXAMPLE_FILES,
    RUN_FILE,
    EOS_TMP,
    GITHUB_ORG,
    S3_BUCKET_URL_ZIP,
    DOCKERHUB_ORG,
)
from ..hub.fetch.actions.template_resolver import TemplateResolver
from ..io.input import ExampleGenerator
from ..hub.content.card import ModelCard
from ..utils.download import GitHubDownloader, S3Downloader
from ..utils.conda import SimpleConda
from ..utils.exceptions_utils import test_exceptions as texc
from ..utils.docker import SimpleDocker
from ..utils.logging import make_temp_dir
from ..utils.spinner import show_loader
from ..utils.hdf5 import Hdf5DataLoader
from ..utils.terminal import run_command_check_output, yes_no_input
from .inspect import ModelInspector
from ..cli import echo


class Options(Enum):
    """
    Enum for different options.
    """

    NUM_SAMPLES = 2
    BASE = "base"
    OUTPUT_CSV = "result.csv"
    INPUT_CSV = "input.csv"
    OUTPUT1_CSV = "output1.csv"
    OUTPUT2_CSV = "output2.csv"
    OUTPUT_FILES = [
        "file.csv",
        "file.json",
        "file.h5",
    ]  # ["file.csv", "file.json", "file.h5"]
    INPUT_TYPES = ["str", "list", "csv"]


class TableType(Enum):
    """
    Enum for different table types.
    """

    MODEL_INFORMATION_CHECKS = "Model Metadata Checks"
    MODEL_FILE_CHECKS = "Model File Checks"
    MODEL_DIRECTORY_SIZES = "Model Directory Sizes"
    MODEL_ENV_SIZES = "Model Environment Sizes"
    RUNNER_CHECKUP_STATUS = "Runner Checkup Status"
    FINAL_RUN_SUMMARY = "Test Run Summary"
    DEPENDECY_CHECK = "Dependency Check"
    COMPUTATIONAL_PERFORMANCE = "Computational Performance Summary"
    SHALLOW_CHECK_SUMMARY = "Validation and Size Check Summary"
    CONSISTENCY_BASH = "Consistency Summary Between Ersilia and Bash Execution Outputs"
    MODEL_OUTPUT = "Model Output Content Validation Summary"
    INSPECT_SUMMARY = "Inspect Summary"


@dataclass
class TableConfig:
    """
    Configuration for a table.
    """

    title: str
    headers: List[str]


TABLE_CONFIGS = {
    TableType.MODEL_INFORMATION_CHECKS: TableConfig(
        title="Model Information Checks", headers=["Check", "Status"]
    ),
    TableType.MODEL_FILE_CHECKS: TableConfig(
        title="Model File Checks", headers=["Check", "Status"]
    ),
    TableType.MODEL_DIRECTORY_SIZES: TableConfig(
        title="Model Directory Sizes", headers=["Check", "Size"]
    ),
    TableType.RUNNER_CHECKUP_STATUS: TableConfig(
        title="Runner Checkup Status",
        headers=["Runner", "Status"],
    ),
    TableType.FINAL_RUN_SUMMARY: TableConfig(
        title="Test Run Summary", headers=["Check", "Status"]
    ),
    TableType.DEPENDECY_CHECK: TableConfig(
        title="Dependency Check", headers=["Check", "Status"]
    ),
    TableType.COMPUTATIONAL_PERFORMANCE: TableConfig(
        title="Computational Performance Summary", headers=["Check", "Status"]
    ),
    TableType.SHALLOW_CHECK_SUMMARY: TableConfig(
        title="Validation and Size Check Results", headers=["Check", "Status"]
    ),
    TableType.MODEL_OUTPUT: TableConfig(
        title="Model Output Content Validation Summary",
        headers=["Check", "Detail", "Status"],
    ),
    TableType.CONSISTENCY_BASH: TableConfig(
        title="Consistency Summary Between Ersilia and Bash Execution Outputs",
        headers=["Check", "Result", "Status"],
    ),
    TableType.INSPECT_SUMMARY: TableConfig(
        title="Inspect Summary", headers=["Check", "Status"]
    ),
}


class STATUS_CONFIGS(Enum):
    """
    Enum for status configurations.
    """

    PASSED = ("PASSED", "green", "✔")
    FAILED = ("FAILED", "red", "✘")
    WARNING = ("WARNING", "yellow", "⚠")
    SUCCESS = ("SUCCESS", "green", "★")
    NA = ("N/A", "dim", "~")

    def __init__(self, label, color, icon):
        self.label = label
        self.color = color
        self.icon = icon

    def __str__(self):
        return f"[{self.color}]{self.icon} {self.label}[/{self.color}]"


# fmt: off
class TestResult(Enum):
    """
    Enum for test results.
    """
    DATE_TIME_RUN = (
        "Date and Time Run",
        lambda: datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    TIME_ELAPSED = (
        "Time to Run Tests (seconds)",
        lambda elapsed: elapsed
    )
    BASIC_CHECKS = (
        "Basic Checks Passed",
        lambda svc: svc.information_check
    )
    SINGLE_INPUT = (
        "Single Input Run Without Error",
        lambda svc: svc.single_input
    )
    EXAMPLE_INPUT = (
        "Example Input Run Without Error",
        lambda svc: svc.example_input
    )
    CONSISTENT_OUTPUT = (
        "Outputs Consistent",
        lambda svc: svc.consistent_output
    )
    BASH_RUN = (
        "Bash Run Without Error",
        lambda run_bash: run_bash,
    )
# fmt: on
    def __init__(self, key, value_function):
        self.key = key
        self.value_function = value_function

    @classmethod
    def generate_results(cls, checkup_service, elapsed_time, run_using_bash):
        """
        Generate test results.

        Parameters
        ----------
        checkup_service : object
            The checkup service.
        elapsed_time : float
            The elapsed time.
        run_using_bash : bool
            Whether to run using bash.

        Returns
        -------
        dict
            The generated results.
        """
        results = {}
        for test in cls:
            func_args = {}
            if "svc" in test.value_function.__code__.co_varnames:
                func_args["svc"] = checkup_service
            if "elapsed" in test.value_function.__code__.co_varnames:
                func_args["elapsed"] = elapsed_time
            if "run_bash" in test.value_function.__code__.co_varnames:
                func_args["run_bash"] = run_using_bash

            value = test.value_function(**func_args)
            results[test.key] = value
        return results

class CheckStrategy:
    """
    Execuetd a strategy for checking inspect commands.

    Parameters
    ----------
    check_function : callable
        The function to check.
    success_key : str
        The key for success.
    details_key : str
        The key for details.
    """
    def __init__(self, check_function, success_key, details_key):
        self.check_function = check_function
        self.success_key = success_key
        self.details_key = details_key

    def execute(self):
        """
        Execute the check strategy.

        Returns
        -------
        dict
            The results of the check.
        """
        if self.check_function is None:
            return {}
        result = self.check_function()
        if result is None:
            return {}
        return {
            self.success_key: result.success,
            self.details_key: result.details,
        }


class InspectService(ErsiliaBase):
    """
    Service for inspecting models and their configurations.

    Parameters
    ----------
    dir : str, optional
        Directory where the model is located.
    model : str, optional
        Model identifier.
    config_json : str, optional
        Path to the configuration JSON file.
    credentials_json : str, optional
        Path to the credentials JSON file.

    Examples
    --------
    .. code-block:: python

        inspector = InspectService(dir="/path/to/model", model="model_id")
        results = inspector.run()
    """

    def __init__(self, dir: str = None, model: str = None, remote: bool = False, config_json: str = None, credentials_json: str = None):
        super().__init__(config_json, credentials_json)
        self.dir = dir

        self.model = model
        self.remote = remote
        self.resolver = TemplateResolver(
            model_id=model,
            repo_path=self.dir
        )

    def run(self, check_keys: list = None) -> dict:
        """
        Run the inspection checks on the specified model.

        Parameters
        ----------
        check_keys : list, optional
            A list of check keys to execute. If None, all checks will be executed.

        Returns
        -------
        dict
            A dictionary containing the results of the inspection checks.

        Raises
        ------
        ValueError
            If the model is not specified.
        KeyError
            If any of the specified keys do not exist.
        """
        def _transform_key(value):
            if value is True:
                return str(STATUS_CONFIGS.PASSED)
            elif value is False:
                return str(STATUS_CONFIGS.FAILED)
            return value

        inspector = ModelInspector(self.model, self.dir)
        checks = self._get_checks(inspector)

        output = {}

        if check_keys:
            for key in check_keys:
                try:
                    strategy = checks.get(key)
                    if strategy.check_function:
                        output.update(strategy.execute())
                except KeyError:
                    raise KeyError(f"Check '{key}' does not exist.")
        else:
            for key, strategy in checks.all():
                if strategy.check_function:
                    output.update(strategy.execute())
        output = {
            " ".join(word.capitalize()
            for word in k.split("_")): _transform_key(v)
            for k, v in output.items()
            }

        output = [(key, value) for key, value in output.items()]
        return output

    def _get_checks(self, inspector: ModelInspector) -> dict:
        def create_check(check_fn, key, details):
            return lambda: CheckStrategy(check_fn, key, details)
        dependency_check = "Dockerfile" if self.resolver.is_bentoml() else "Install_YAML"
        checks = {
            "is_github_url_available": create_check(
                inspector.check_repo_exists if self.remote else lambda: None,
                "is_github_url_available",
                "is_github_url_available_details",
            ),
            "complete_metadata": create_check(
                inspector.check_complete_metadata if self.remote else lambda: None,
                "complete_metadata",
                "complete_metadata_details",
            ),
            "complete_folder_structure": create_check(
                inspector.check_complete_folder_structure,
                "complete_folder_structure",
                "complete_folder_structure_details",
            ),
            "docker_check": create_check(
                inspector.check_dependencies_are_valid,
                f"{dependency_check}_check",
                "check_details",
            ),
            "computational_performance_tracking": create_check(
                inspector.check_computational_performance,
                "computational_performance_tracking",
                "computational_performance_tracking_details",
            ),
            "extra_files_check": create_check(
                inspector.check_no_extra_files if self.remote else lambda: None,
                "extra_files_check",
                "extra_files_check_details",
            ),
        }

        class LazyChecks:
            def __init__(self, checks):
                self._checks = checks
                self._loaded = {}

            def get(self, key):
                if key not in self._loaded:
                    if key not in self._checks:
                        raise KeyError(f"Check '{key}' does not exist.")
                    self._loaded[key] = self._checks[key]()
                return self._loaded[key]

            def all(self):
                for key in self._checks.keys():
                    yield key, self.get(key)

        return LazyChecks(checks)


class SetupService:
    """
    Service for setting up the environment and fetching the model repository.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    dir : str
        Directory where the model repository will be cloned.
    logger : logging.Logger
        Logger for logging messages.
    from_github : bool
        Flag indicating whether to fetch the repository from GitHub.
    from_s3 : bool
        Flag indicating whether to fetch the repository from S3.
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
        self.s3 = (
            self.metadata.get("card", {}).get("S3") or
            self.metadata.get("metadata", {}).get("S3")
        )
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
        if not self.overwrite:
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
            zip_ref.extractall(self.dir)

    def _download_github(self):
        if not self._handle_overwrite():
            self.logger.info("Skipping GitHub download as user chose not to overwrite.")
            return

        destination = os.path.join(self.dir, self.model_id)
        self.logger.info(f"Cloning repository from GitHub to: {destination}")
        self.github_down.clone(
            org=GITHUB_ORG,
            repo=self.model_id,
            destination=destination,
        )

    def get_model(self):
        if self.from_s3:
            self._download_s3()

        if self.from_github:
            self._download_github()
                
    @staticmethod
    def run_command(command: str, logger, capture_output: bool = False, shell: bool = True) -> str:
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
                    shell=shell
                )
                return result.stdout
            else:
                process = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    shell=shell
                )

                stdout_lines, stderr_lines = [], []

                for line in iter(process.stdout.readline, ''):
                    if line.strip():
                        stdout_lines.append(line.strip())
                        logger.info(line.strip())

                for line in iter(process.stderr.readline, ''):
                    if line.strip():
                        stderr_lines.append(line.strip())
                        logger.info(line.strip())

                process.wait()
                if process.returncode != 0:
                    raise subprocess.CalledProcessError(
                        returncode=process.returncode,
                        cmd=command,
                        output="\n".join(stdout_lines),
                        stderr="\n".join(stderr_lines),
                    )

                return process.stdout

        except subprocess.CalledProcessError as e:
            logger.error(f"Error executing command: {e}")
            if e.output:
                logger.debug(f"Output: {e.output.strip()}")
            if e.stderr:
                logger.debug(f"Error: {e.stderr.strip()}")
            sys.exit(1)
        except Exception as e:
            logger.debug(f"Unexpected error: {e}")
            sys.exit(1)

    def check_conda_env(self):
        """
        Check if the Conda environment for the model exists.

        Raises
        ------
        Exception
            If the Conda environment does not exist.
        """
        if self.conda.exists(self.model_id):
            self.logger.debug(
                f"Conda environment '{self.model_id}' already exists."
            )
        else:
            raise Exception(
                f"Conda virtual environment not found for {self.model_id}"
            )

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
                "conda env list",
                logger=logger,
                capture_output=True
            )
            for line in result.splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if parts[0] == model_id:
                    return parts[-1]
        except subprocess.CalledProcessError as e:
            print(f"Error running conda command: {e.stderr}")
        except Exception as e:
            print(f"Unexpected error: {e}")

        return None

class IOService:
    """
    Service for handling input/output operations related to model testing.

    Parameters
    ----------
    logger : logging.Logger
        Logger for logging messages.
    dest_dir : str
        Destination directory for storing model-related files.
    model_path : str
        Path to the model.
    bundle_path : str
        Path to the model bundle.
    bentoml_path : str
        Path to the BentoML files.
    model_id : str
        Identifier of the model.
    dir : str
        Directory where the model repository is located.

    Examples
    --------
    .. code-block:: python

        ios = IOService(logger=logger, dest_dir="/path/to/dest", model_path="/path/to/model",
                        bundle_path="/path/to/bundle", bentoml_path="/path/to/bentoml",
                        model_id="model_id", dir="/path/to/dir")
        ios.read_information()
    """

    # Required files
    RUN_FILE = f"model/framework/{RUN_FILE}"
    BENTOML_FILES = [
        DOCKERFILE_FILE,
        METADATA_JSON_FILE,
        RUN_FILE,
        "src/service.py",
        "pack.py",
        "README.md",
        "LICENSE"
    ]

    ERSILIAPACK_FILES = [
        INSTALL_YAML_FILE,
        METADATA_YAML_FILE,
        PREDEFINED_EXAMPLE_FILES[0],
        PREDEFINED_EXAMPLE_FILES[1],
        RUN_FILE,
        "README.md",
        "LICENSE",
    ]

    def __init__(self, logger,  model_id: str, dir: str):
        self.logger = logger
        self.model_id = model_id
        self.dir = dir
        self.model_size = 0
        self.console = Console()
        self.check_results = []
        self.simple_docker = SimpleDocker()

    def _run_check(self, check_function, data, check_name: str, additional_info=None) -> bool:
        try:
            if additional_info is not None:
                check_function(additional_info)
            else:
                check_function(data)
            self.check_results.append((
                check_name,
                str(STATUS_CONFIGS.PASSED)
            ))
            return True
        except Exception as e:
            self.logger.error(
                f"Check '{check_name}' failed: {e}"
            )
            self.check_results.append((
                check_name,
                str(STATUS_CONFIGS.FAILED)
            ))
            return False

    def _generate_table(self, title: str, headers: List[str], rows: List[List[str]], large_table: bool = False, merge: bool = False):
        f_col_width = 30 if large_table else 30
        l_col_width = 50 if large_table else 10
        d_col_width = 30 if not large_table else 20

        table = Table(
            title=Text(
                title,
                style="bold light_green"
            ),
            border_style="light_green",
            show_lines=True,
        )

        table.add_column(
            headers[0],
            justify="left",
            width=f_col_width,
            style="bold"
        )
        for header in headers[1:-1]:
            table.add_column(
                header,
                justify="center",
                width=d_col_width,
                style="bold"
            )
        table.add_column(
            headers[-1],
            justify="right",
            width=l_col_width,
            style="bold"
        )

        prev_value = None
        for row in rows:
            first_col = str(row[0])
            if merge and first_col == prev_value:
                first_col = ""
            else:
                prev_value = first_col

            styled_row = [
                Text(first_col, style="bold"),
                *[str(cell) for cell in row[1:-1]],
                row[-1]
            ]
            table.add_row(*styled_row)

        self.console.print(table)

    @staticmethod
    def get_model_type(model_id: str, repo_path: str) -> str:
        """
        Get the type of the model based on the repository contents.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        repo_path : str
            Path to the model repository.

        Returns
        -------
        str
            The type of the model (e.g., PACK_METHOD_BENTOML, PACK_METHOD_FASTAPI).
        """
        resolver = TemplateResolver(
            model_id=model_id,
            repo_path=repo_path
        )
        if resolver.is_bentoml():
            return PACK_METHOD_BENTOML
        elif resolver.is_fastapi():
            return PACK_METHOD_FASTAPI
        else:
            return None

    def get_file_requirements(self) -> List[str]:
        """
        Get the list of required files for the model.

        Returns
        -------
        List[str]
            List of required files.

        Raises
        ------
        ValueError
            If the model type is unsupported.
        """
        type = IOService.get_model_type(
            model_id=self.model_id,
            repo_path=self.dir
        )
        if type == PACK_METHOD_BENTOML:
            return self.BENTOML_FILES
        elif type == PACK_METHOD_FASTAPI:
            return self.ERSILIAPACK_FILES
        else:
            raise ValueError(
                f"Unsupported model type: {type}"
            )

    def read_information(self) -> dict:
        """
        Read the information file for the model.

        Returns
        -------
        dict
            The contents of the information file.

        Raises
        ------
        FileNotFoundError
            If the information file does not exist.
        """
        file = os.path.join(
            EOS_TMP,
            self.model_id,
            METADATA_JSON_FILE
        )
        if not os.path.exists(file):
            raise FileNotFoundError(
                f"Information file does not exist for model {self.model_id}"
        )
        with open(file, "r") as f:
            return json.load(f)

    def print_output(self, result, output):
        """
        Print the output of a result.

        Parameters
        ----------
        result : any
            The result to print.
        output : file-like object
            The output file to write to.
        """
        def write_output(data):
            if output is not None:
                with open(output.name, "w") as file:
                    json.dump(data, file)
            else:
                self.logger.debug(json.dumps(data, indent=4))

        if isinstance(result, types.GeneratorType):
            for r in result:
                write_output(r if r is not None else "Something went wrong")
        else:
            self.logger.debug(result)

    def get_conda_env_size(self) -> int:
        """
        Get the size of the Conda environment for the model.

        Returns
        -------
        int
            The size of the Conda environment in megabytes.

        Raises
        ------
        Exception
            If there is an error calculating the size.
        """
        try:
            loc = SetupService.get_conda_env_location(
                self.model_id,
                self.logger
            )
            return self.calculate_directory_size(loc)
        except Exception as e:
            self.logger.error(
                f"Error calculating size of Conda environment '{self.model_id}': {e}"
            )
            return 0

    def calculate_directory_size(self, path: str) -> int:
        """
        Calculate the size of a directory.

        Parameters
        ----------
        path : str
            The path to the directory.

        Returns
        -------
        int
            The size of the directory.
        """
        try:
            size_output = SetupService.run_command(
                ["du", "-sm", path],
                logger=self.logger,
                capture_output=True,
                shell=False
            )
            size = int(size_output.split()[0])
            return size
        except Exception as e:
            self.logger.error(
                f"Error calculating directory size for {path}: {e}"
            )
            return 0
        
    def calculate_image_size(self):
        """
        Calculates the size of a Docker image.

        Args:
            image_name (str): The name (and optionally tag) of the Docker image, e.g., "ubuntu:latest".

        Returns:
            str: The size of the image in a human-readable format, or an error message if the image is not found.
        """
        image_name = f"{DOCKERHUB_ORG}/{self.model_id}"
        client = docker.from_env()  # Initialize Docker client
        try:
            image = client.images.get(image_name)
            size_bytes = image.attrs['Size']  # Get the size in bytes
            size_mb = size_bytes / (1024 ** 2)  # Convert bytes to MB
            return f"{size_mb:.2f} MB"
        except docker.errors.ImageNotFound:
            return f"Image '{image_name}' not found."
        except Exception as e:
            return f"An error occurred: {e}"



    @throw_ersilia_exception()
    def get_directories_sizes(self) -> str:
        """
        Get the sizes of the model directory.

        Returns
        -------
        str
          A string of containing size of the model directory
        """
        dir_size = self.calculate_directory_size(self.dir)
        dir_size = f"{dir_size}MB"
        return dir_size

    @throw_ersilia_exception()
    def get_env_sizes(self) -> str:
        """
        Get the sizes of the Conda environment directory.

        Returns
        -------
        str
          A string of containing size of the model environment
        """
        env_size = self.get_conda_env_size()
        env_size = f"{env_size}MB"
        return env_size

class CheckService:
    """
    Service for performing various checks on the model.

    Parameters
    ----------
    logger : logging.Logger
        Logger for logging messages.
    model_id : str
        Identifier of the model.
    dest_dir : str
        Destination directory for storing model-related files.
    dir : str
        Directory where the model repository is located.
    ios : IOService
        Instance of IOService for handling input/output operations.

    Examples
    --------
    .. code-block:: python

        check_service = CheckService(logger=logger, model_id="model_id", dest_dir="/path/to/dest",
                                     dir="/path/to/dir", ios=ios)
        check_service.check_files()
    """

    # model specs
    MODEL_TASKS = {
        "Classification",
        "Regression",
        "Generative",
        "Representation",
        "Similarity",
        "Clustering",
        "Dimensionality reduction",
    }

    MODEL_OUTPUT = {
        "Boolean",
        "Compound",
        "Descriptor",
        "Distance",
        "Experimental value",
        "Image",
        "Other value",
        "Probability",
        "Protein",
        "Score",
        "Text",
    }

    INPUT_SHAPE = {
        "Single",
        "Pair",
        "List",
        "Pair of Lists",
        "List of Lists"
    }

    OUTPUT_SHAPE = {
        "Single",
        "List",
        "Flexible List",
        "Matrix",
        "Serializable Object"
    }

    def __init__(self, logger, model_id: str, dir: str, ios: IOService):
        self.logger = logger
        self.model_id = model_id
        self.dir = dir
        self._run_check = ios._run_check
        self._generate_table = ios._generate_table
        self._print_output = ios.print_output
        self.get_file_requirements = ios.get_file_requirements
        self.console = ios.console
        self.check_results = ios.check_results
        self.resolver = TemplateResolver(
            model_id=model_id,
            repo_path=self.dir
        )
    
    def _get_metadata(self):
        path = METADATA_JSON_FILE if self.resolver.is_bentoml() else METADATA_YAML_FILE
        path = os.path.join(self.dir, path)

        with open(path, 'r') as file:
            if path.endswith('.json'):
                data = json.load(file)
            elif path.endswith(('.yml', '.yaml')):
                data = yaml.safe_load(file)
            else:
                raise ValueError(f"Unsupported file format: {path}")
        return data

    def _check_file_existence(self, path):
        if not os.path.exists(os.path.join(self.dir, path)):
            raise FileNotFoundError(
                f"File '{path}' does not exist."
        )

    def check_files(self):
        """
        Check the existence of required files for the model.
        """
        requirements = self.get_file_requirements()
        for file in requirements:
            self.logger.debug(f"Checking file: {file}")
            self._run_check(
                self._check_file_existence,
                None,
                f"File: {file}",
                file
            )

    def _check_model_id(self, data):
        self.logger.debug("Checking model ID...")
        if data["Identifier"] != self.model_id:
            raise texc.WrongCardIdentifierError(self.model_id)


    def _check_model_slug(self, data):
        self.logger.debug("Checking model slug...")
        if not data["Slug"]:
            raise texc.EmptyField("slug")


    def _check_model_description(self, data):
        self.logger.debug("Checking model description...")
        if not data["Description"]:
            raise texc.EmptyField("Description")
    
    def _check_model_tag(self, data):
        self.logger.debug("Checking model tag...")
        if not data["Tag"]:
            raise texc.EmptyField("Tag")
        
    def _check_model_source_code(self, data):
        self.logger.debug("Checking model source code...")
        if not data["Source Code"]:
            raise texc.EmptyField("Source Code")
        
    def _check_model_source_title(self, data):
        self.logger.debug("Checking model title...")
        if not data["Title"]:
            raise texc.EmptyField("Title")
        
    def _check_model_status(self, data):
        self.logger.debug("Checking model status...")
        if not data["Status"]:
            raise texc.EmptyField("Status")
        
    def _check_model_contributor(self, data):
        self.logger.debug("Checking model contributor...")
        if not data["Contributor"]:
            raise texc.EmptyField("Contributor") 
        
    def _check_model_interpret(self, data):
        self.logger.debug("Checking model interpretation...")
        if not data["Interpretation"]:
            raise texc.EmptyField("Interpretation")
        
    def _check_model_dockerhub_url(self, data):
        key = "DockerHub"
        self.logger.info(f"Data: {data}")   
        self.logger.debug(f"Checking {key} URL field..")
        if key in data:
            self.logger.debug(f"Checking {key} URL field..")
            if not data[key]:
                self.logger.debug(f"Checking {key} URL field..")
                raise texc.EmptyField(key)
        else:
            self.logger.debug(f"Checking {key} URL field..")
            raise texc.EmptyKey(key)
        
    def _check_model_s3_url(self, data):
        key = "S3"
        self.logger.debug(f"Checking {key} URL field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)
        
    def _check_model_arch(self, data):
        key = "Docker Architecture"
        self.logger.debug(f"Checking {key} field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)
        
    def _check_model_publication(self, data):
        key = "Publication"
        self.logger.debug(f"Checking {key} field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)
        
    def _check_model_task(self, data):
        self.logger.debug("Checking model task...")
        raw_tasks = data.get("Task", "")
        if isinstance(raw_tasks, str):
            tasks = [
                task.strip()
                for task
                in raw_tasks.split(",")
                if task.strip()
            ]
        elif isinstance(raw_tasks, list):
            tasks = [
                task.strip()
                for task
                in raw_tasks
                if isinstance(task, str) and task.strip()
            ]
        else:
            raise texc.InvalidEntry(
                "Task",
                message="Task field must be a string or list."
            )

        if not tasks:
            raise texc.InvalidEntry(
                "Task",
                message="Task field is missing or empty."
            )

        invalid_tasks = [task for task in tasks if task not in self.MODEL_TASKS]
        if invalid_tasks:
            raise texc.InvalidEntry(
                "Task", message=f"Invalid tasks: {', '.join(invalid_tasks)}"
            )

        self.logger.debug("All tasks are valid.")

    def _check_model_output(self, data):
        self.logger.debug("Checking model output...")
        raw_outputs = data.get("Output", "")
        if isinstance(raw_outputs, str):
            outputs = [
                output.strip()
                for output
                in raw_outputs.split(",")
                if output.strip()
            ]
        elif isinstance(raw_outputs, list):
            outputs = [
                output.strip()
                for output
                in raw_outputs
                if isinstance(output, str) and output.strip()
            ]
        else:
            raise texc.InvalidEntry(
                "Output",
                message="Output field must be a string or list."
            )

        if not outputs:
            raise texc.InvalidEntry(
                "Output",
                message="Output field is missing or empty."
            )

        invalid_outputs = [
            output
            for output
            in outputs
            if output not in self.MODEL_OUTPUT
        ]
        if invalid_outputs:
            raise texc.InvalidEntry(
                "Output",
                message=f"Invalid outputs: {' '.join(invalid_outputs)}"
            )

        self.logger.debug("All outputs are valid.")

    def _check_model_input(self, data):
        self.logger.debug("Checking model input")
        valid_inputs = [{"Compound"}, {"Protein"}, {"Text"}]

        model_input = data.get("Input")

        if not model_input or set(model_input) not in valid_inputs:
            raise texc.InvalidEntry("Input")

    def _check_model_input_shape(self, data):
        self.logger.debug("Checking model input shape")
        model_input_shape = data.get("Input Shape")

        if model_input_shape not in self.INPUT_SHAPE:
            raise texc.InvalidEntry("Input Shape")

    def _check_model_output_type(self, data):
        self.logger.debug("Checking model output type...")
        valid_output_types = [{"String"}, {"Float"}, {"Integer"}]

        model_output_type = (
            data.get("Output Type")
        )

        if not model_output_type or set(model_output_type) not in valid_output_types:
            raise texc.InvalidEntry("Output Type")

    def _check_model_output_shape(self, data):
        self.logger.debug("Checking model output shape...")
        model_output_shape = (
            data.get("Output Shape")
        )

        if model_output_shape not in self.OUTPUT_SHAPE:
            raise texc.InvalidEntry("Output Shape")

    @throw_ersilia_exception()
    def check_information(self):
        """
        Perform various checks on the model information.

        Parameters
        ----------
        output : file-like object
            The output file to write to.
        """
        self.logger.debug(f"Beginning checks for {self.model_id} model information")
        data = self._get_metadata()

        self._run_check(self._check_model_id, data, "Model ID")
        self._run_check(self._check_model_slug, data, "Model Slug")
        self._run_check(self._check_model_status, data, "Model Status")
        self._run_check(self._check_model_source_title, data, "Model Title")
        self._run_check(self._check_model_description, data, "Model Description")
        self._run_check(self._check_model_task, data, "Model Task")
        self._run_check(self._check_model_input, data, "Model Input")
        self._run_check(self._check_model_input_shape, data, "Model Input Shape")
        self._run_check(self._check_model_output, data, "Model Output")
        self._run_check(self._check_model_output_type, data, "Model Output Type")
        self._run_check(self._check_model_output_shape, data, "Model Output Shape")
        self._run_check(self._check_model_interpret, data, "Model Interpretation")
        self._run_check(self._check_model_tag, data, "Model Tag")
        self._run_check(self._check_model_publication, data, "Model Publication")
        self._run_check(self._check_model_source_code, data, "Model Source Code")
        self._run_check(self._check_model_contributor, data, "Model Contributor")
        self._run_check(self._check_model_dockerhub_url, data, "Model Dockerhub URL")
        self._run_check(self._check_model_s3_url, data, "Model S3 URL")
        self._run_check(self._check_model_arch, data, "Model Docker Architecture")
    
    def get_inputs(self, run_example, types):
        samples = run_example(
            n_samples=Options.NUM_SAMPLES.value,
            file_name=None,
            simple=True,
            try_predefined=False
        )
        samples = [sample["input"] for sample in samples]
        if types == "str":
            return samples[0]
        if types == "list":
            return json.dumps(samples)
        if types == "csv":
            run_example(
            n_samples=Options.NUM_SAMPLES.value,
            file_name=Options.INPUT_CSV.value,
            simple=True,
            try_predefined=False
            )
            return Options.INPUT_CSV.value

    def check_model_output_content(self, run_example, run_model):
        status = []
        self.logger.debug("Checking model output...")
        for inp_type in Options.INPUT_TYPES.value:
            for i, output_file in enumerate(Options.OUTPUT_FILES.value):
                self.logger.debug(f"Checking model output for input type: {inp_type} and output file: {output_file}")
                inp_data = self.get_inputs(run_example, inp_type)
                self.logger.debug(f"Input data: {inp_data}")
                run_model(inputs=inp_data, output=output_file, batch=100)
                _status = self.validate_file_content(output_file, inp_type)
                status.append(_status)
                self.logger.debug(f"Output file {output_file} is {_status}")
        return status
    
    def validate_file_content(self, file_path, input_type):
        def is_invalid_value(item):
            return item in [None, "", "null"] or (isinstance(item, float) and math.isnan(item))
        
        def validate_json(content):
            data_structures = {
                "Single": lambda x: isinstance(x, list) and len(x) == 1 or isinstance(x, (int, float, str)),
                "List": lambda x: isinstance(x, list) and len(x) > 1 and all(isinstance(item, (int, float, str)) for item in x),
                "Flexible List": lambda x: isinstance(x, list) and all(isinstance(item, (str, int, float)) for item in x),
                "Matrix": lambda x: isinstance(x, list) and all(
                    isinstance(row, list) and all(isinstance(item, (int, float)) for item in row) for row in x
                ),
                "Serializable Object": lambda x: isinstance(x, dict),
            }
            
            for check_type, validator in data_structures.items():
                types = validator(content)
                self.logger.debug(f"Checking outcome of type '{types}'...")
                if validator(content):
                    content = content[0] if isinstance(content, list) and len(content) == 1 else content
                    if check_type == "Single" and is_invalid_value(content):
                        self.logger.debug(f"Invalid value '{content}' in outcome of type 'Single'.")
                        raise ValueError(f"Invalid value '{content}' in outcome of type 'Single'.")
                    elif check_type in ["List", "Flexible List"]:
                        self.logger.debug(f"Checking outcome of type '{check_type}'...")
                        for item in content:
                            if is_invalid_value(item):
                                raise ValueError(f"Invalid value '{item}' in outcome of type '{check_type}'.")
                    elif check_type == "Matrix":
                        self.logger.debug(f"Checking outcome of type 'Matrix'...")
                        for row in content:
                            for item in row:
                                if is_invalid_value(item):
                                    raise ValueError(f"Invalid value '{item}' in outcome of type 'Matrix'.")
                    elif check_type == "Serializable Object":
                        self.logger.debug(f"Checking outcome of type 'Serializable Object'...")
                        for key, value in content.items():
                            if is_invalid_value(value):
                                raise ValueError(f"Invalid value '{value}' for key '{key}' in outcome of type 'Serializable Object'.")
                    return f"{input_type} : JSON", "Valid Content", str(STATUS_CONFIGS.PASSED)
        
            raise TypeError("Unknown content structure.")
        
        def check_json(file_path):
            with open(file_path, "r") as f:
                try:
                    content = json.load(f)
                    self.logger.debug(f"Checking JSON file: {file_path}")
                    self.logger.debug(f"Content: {content}")
                    for item in content:
                        # outcome = item.get("output", {}).get("value", None)
                        outcome = next(iter(item.values()), None)
                        self.logger.debug(f"Checking outcome: {outcome}")
                        if outcome:
                            return validate_json(outcome)
                    return f"{input_type} : Unknown", "No valid JSON outcome found", str(STATUS_CONFIGS.FAILED)
                except json.JSONDecodeError as e:
                    return f"{input_type} : JSON", f"Invalid JSON content: {e}", str(STATUS_CONFIGS.FAILED)
        
        def check_csv(file_path):
            self.logger.debug(f"Checking CSV file: {file_path}")
            with open(file_path, "r") as f:
                reader = csv.reader(f)
                rows = list(reader)[1:] 
                if any(any(is_invalid_value(cell) for cell in row) for row in rows):
                    return f"{input_type} : CSV", "Invalid values found in CSV content", str(STATUS_CONFIGS.FAILED)
                return f"{input_type} : CSV", "Valid Content", str(STATUS_CONFIGS.PASSED)
            
        def check_h5(file_path):
            self.logger.debug(f"Checking HDF5 file: {file_path}")
            try:
                loader = Hdf5DataLoader()
                loader.load(file_path)
                content = next((x for x in [loader.values, loader.keys, loader.inputs, loader.features] if x is not None), None)
                
                self.logger.debug(f"Content: {content}")    
                if content is None or (hasattr(content, 'size') and content.size == 0):
                    return f"{input_type} : HDF5", "Empty content", str(STATUS_CONFIGS.FAILED)
                
                if any(any(is_invalid_value(cell) for cell in row) for row in content):
                    return f"{input_type} : HDF5", "Invalid values found in HDF5 content"
                
                return f"{input_type} : HDF5", "Valid Content", str(STATUS_CONFIGS.PASSED)
            except Exception as e:
                return f"{input_type} : HDF5", f"Invalid HDF5 content: {e}", str(STATUS_CONFIGS.FAILED)

        
        if not Path(file_path).exists():
            raise FileNotFoundError(f"File {file_path} does not exist.")
        
        file_extension = Path(file_path).suffix.lower()
        if file_extension == ".json":
            return check_json(file_path)
        elif file_extension == ".csv":
            return check_csv(file_path)
        elif file_extension == ".h5":
            return check_h5(file_path)
        else:
            raise ValueError(f"Unsupported file type: {file_extension}. Supported types are JSON, CSV, and HDF5.")
        
    @throw_ersilia_exception()
    def check_single_input(self, run_model, run_example):
        """
        Check if the model can run with a single input to check if it has a value
        in the produced output csv.

        Parameters
        ----------
        output : file-like object
            The output file to write to.
        run_model : callable
            Function to run the model.
        run_example : callable
            Function to generate example input.
        """
        self.logger.debug("Checking model with single input...")
        output = Options.OUTPUT_CSV.value
        input = run_example(
            n_samples=Options.NUM_SAMPLES.value,
            file_name=None,
            simple=True,
            try_predefined=False
        )
        input = json.dumps([input['input'] for input in input])
        run_model(
            inputs=input,
            output=output,
            batch=100
        )

        def read_csv(file_path):
            absolute_path = os.path.abspath(file_path)
            if not os.path.exists(absolute_path):
                raise FileNotFoundError(f"File not found: {absolute_path}")
            with open(absolute_path, mode='r') as csv_file:
                reader = csv.DictReader(csv_file)
                return [row for row in reader]

        csv_content, _completed_status = read_csv(output), []
        if csv_content:
            _completed_status.append(("Check Single Input", str(STATUS_CONFIGS.PASSED)))
        else:
            _completed_status.append(("Check Single Input", str(STATUS_CONFIGS.FAILED)))
        return _completed_status
        
    @throw_ersilia_exception()
    def check_example_input(self, run_model, run_example):
        """
        Check if the model can run with example input without error.

        Parameters
        ----------
        output : file-like object
            The output file to write to.
        run_model : callable
            Function to run the model.
        run_example : callable
            Function to generate example input.
        """
        self.logger.debug("Checking model with example input...")
        output = Options.OUTPUT_CSV.value
        input = run_example(
            n_samples=Options.NUM_SAMPLES.value,
            file_name=None,
            simple=True,
            try_predefined=True
        )
        input = json.dumps([input['input'] for input in input])

        self.logger.debug("Testing model on input of 5 smiles given by 'example' command")

        run_model(
            inputs=input,
            output=output,
            batch=100
        )

        csv_content, _completed_status = input, []
        if csv_content:
            _completed_status.append(("Check Predefined Example Input", str(STATUS_CONFIGS.PASSED)))
        else:
            _completed_status.append(("Check Predefined Example Input", str(STATUS_CONFIGS.FAILED)))
        return _completed_status

    @throw_ersilia_exception()
    def check_consistent_output(self, run_example, run_model):
        """
        Check if the model produces consistent output.

        Parameters
        ----------
        run_example : callable
            Function to generate example input.
        run_model : callable
            Function to run the model.
        """
        self.logger.debug("Confirming model produces consistent output...")
        def compute_rmse(y_true, y_pred):
            return sum((yt - yp) ** 2 for yt, yp in zip(y_true, y_pred)) ** 0.5 / len(y_true)

        def _compare_output_strings(output1, output2):
            if output1 is None and output2 is None:
                return 100
            else:
                return fuzz.ratio(output1, output2)

        def validate_output(output1, output2):
            if not isinstance(output1, type(output2)):
                raise texc.InconsistentOutputTypes(self.model_id)

            if output1 is None:
                return

            if isinstance(output1, (float, int)):
                rmse = compute_rmse([output1], [output2])
                if rmse > 0.1:
                    raise texc.InconsistentOutputs(self.model_id)

                rho, _ = spearmanr([output1], [output2])
                if rho < 0.5:
                    raise texc.InconsistentOutputs(self.model_id)

            elif isinstance(output1, list):
                rmse = compute_rmse(output1, output2)
                if rmse > 0.1:
                    raise texc.InconsistentOutputs(self.model_id)

                rho, _ = spearmanr(output1, output2)
                if rho < 0.5:
                    raise texc.InconsistentOutputs(self.model_id)

            elif isinstance(output1, str):
                if _compare_output_strings(output1, output2) <= 95:
                    raise texc.InconsistentOutputs(self.model_id)

        def read_csv(file_path):
            absolute_path = os.path.abspath(file_path)
            if not os.path.exists(absolute_path):
                raise FileNotFoundError(f"File not found: {absolute_path}")
            with open(absolute_path, mode='r') as csv_file:
                reader = csv.DictReader(csv_file)
                return [row for row in reader]

        output1_path = os.path.abspath(Options.OUTPUT1_CSV.value)
        output2_path = os.path.abspath(Options.OUTPUT2_CSV.value)

        self.logger.debug("Confirming model produces consistent output...")

        input = run_example(
            n_samples=Options.NUM_SAMPLES.value,
            file_name=None,
            simple=True,
            try_predefined=False
        )
        input = json.dumps([input['input'] for input in input])

        run_model(
            inputs=input,
            output=output1_path,
            batch=100
        )
        run_model(
            inputs=input,
            output=output2_path,
            batch=100
        )

        data1 = read_csv(output1_path)
        data2 = read_csv(output1_path)
        _completed_status = []
        try:
            for res1, res2 in zip(data1, data2):
                for key in res1:
                    if key in res2:
                        validate_output(res1[key], res2[key])
                    else:
                        raise KeyError(f"Key '{key}' not found in second result.")
            _completed_status.append(("Check Consistency of Model Output", str(STATUS_CONFIGS.PASSED)))
        except:
            _completed_status.append(("Check Consistency of Model Output", f"{str(STATUS_CONFIGS.FAILED)}: Inconsistent Output Detected!"))
            raise texc.InconsistentOutputs("Inconsistent Output Detected!")
        return _completed_status

class RunnerService:
    """
    Service for running model tests and checks.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    logger : logging.Logger
        Logger for logging messages.
    ios_service : IOService
        Instance of IOService for handling input/output operations.
    checkup_service : CheckService
        Instance of CheckService for performing various checks on the model.
    setup_service : SetupService
        Instance of SetupService for setting up the environment and fetching the model repository.
    model_path : str
        Path to the model.
    level : str
        Level of checks to perform.
    dir : str
        Directory where the model repository is located.
    remote : bool
        Flag indicating whether to fetch the repository from a remote source.
    inspect : bool
        Flag indicating whether to perform inspection checks.
    remove : bool
        Flag indicating whether to remove the model directory after tests.
    inspecter : InspectService
        Instance of InspectService for inspecting models and their configurations.
    """

    def __init__(
        self,
        model_id: str,
        logger,
        ios_service: IOService,
        checkup_service: CheckService,
        setup_service: SetupService,
        level: str,
        dir: str,
        model_path: Callable,
        from_github: bool,
        from_s3: bool,
        from_dockerhub: bool,
        version: str,
        shallow: bool,
        deep: bool,
        inspecter: InspectService
    ):
        self.model_id = model_id
        self.logger = logger
        self.setup_service = setup_service
        self.ios_service = ios_service
        self.console = ios_service.console
        self.checkup_service = checkup_service
        self.model_path = model_path(self.model_id)
        self.level = level
        self.dir = dir
        self.from_github = from_github
        self.from_s3 = from_s3
        self.from_dockerhub = from_dockerhub
        self.version = version
        self.shallow = shallow
        self.deep = deep
        self.inspecter = inspecter
        self.example = ExampleGenerator(
            model_id=self.model_id
        )
        self.run_using_bash = False

    def run_model(self, inputs: list, output: str, batch: int):
        """
        Run the model with the given input and output parameters.

        Parameters
        ----------
        input : list
            List of input samples.
        output : str
            Path to the output file.
        batch : int
            Batch size for running the model.

        Returns
        -------
        str
            The output of the command.
        """
        self.logger.info(f"Input after escape: {inputs}")
        out = SetupService.run_command(
            f"ersilia -v serve {self.model_id} && ersilia  -v run -i '{inputs}' -o {output} -b {str(batch)}",
            logger=self.logger,
        )
        self.logger.info(out)
        return out
    
    def fetch(self):
        """
        Fetch the model repository from the specified directory.
        """
        def _fetch(dir, model_id, logger):
            loc = (
                ["--from_dir", self.dir]
                if self.from_github or self.from_s3
                else ["--from_dockerhub"] + (["--version", self.version] if self.version else [])
            )
            self.logger.info(f"Fetching model from: {loc}")
            out =SetupService.run_command(
                " ".join(["ersilia",
                "-v",
                "fetch", model_id,
                *loc
                ]),
                logger=logger,
            )
            logger.info(f"Fetch out: {out}")

        if os.path.exists(self.model_path):
            SetupService.run_command(
                " ".join(["ersilia",
                "-v",
                "delete", self.model_id,
                ]),
                logger=self.logger,

            )
        _fetch(self.dir, self.model_id, self.logger)
        
    def run_exampe(
        self,
        n_samples: int,
        file_name: str = None,
        simple: bool = True,
        try_predefined: bool = False
    ):
        """
        Generate example input samples for the model.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate.
        file_name : str, optional
            Name of the file to save the samples.
        simple : bool, optional
            Flag indicating whether to generate simple samples.
        try_predefined : bool, optional
            Flag indicating whether to try predefined samples.

        Returns
        -------
        list
            List of generated input samples.
        """
        return self.example.example(
            n_samples=n_samples,
            file_name=file_name,
            simple=simple,
            try_predefined=try_predefined
        )
    
    @throw_ersilia_exception()
    def run_bash(self):
        """
        Run the model using a bash script and compare the outputs for consistency.

        Raises
        ------
        RuntimeError
            If there is an error during the subprocess execution or output comparison.
        """
        def compute_rmse(y_true, y_pred):
            return sum((yt - yp) ** 2 for yt, yp in zip(y_true, y_pred)) ** 0.5 / len(y_true)
        
        def compare_outputs(bsh_data, ers_data):
            _completed_status = []
            self.logger.debug(f"Bash Data: {bsh_data}")
            self.logger.debug(f"Ers Data: {ers_data}")
            columns = set(bsh_data[0].keys()) & set(data[0].keys())
            self.logger.debug(f"Common columns: {columns}")

            for column in columns:
                bv = [row[column] for row in bsh_data]
                ev = [row[column] for row in ers_data]

                if all(isinstance(val, (int, float)) for val in bv + ev):
                    rmse = compute_rmse(bv, ev)
                    self.logger.debug(f"RMSE for {column}: {rmse}")

                    if rmse > 0.1:
                        _completed_status.append(("RMSE", f"RMSE > 10% | {rmse * 100}%", str(STATUS_CONFIGS.FAILED)))
                        raise texc.InconsistentOutputs(self.model_id)
                    _completed_status.append(("RMSE", f"{rmse * 100}%", str(STATUS_CONFIGS.PASSED)))
                    
                elif all(isinstance(val, str) for val in bv + ev):
                    if not all(
                        self._compare_string_similarity(a, b, 95)
                        for a, b
                        in zip(bv, ev)
                    ):
                        _completed_status.append(("String Similarity", "< 95%", str(STATUS_CONFIGS.FAILED)))
                        raise texc.InconsistentOutputs(self.model_id)
                    _completed_status.append(("String Similarity", "> 95%", str(STATUS_CONFIGS.PASSED)))

            return _completed_status
                    

        def read_csv(path, flag=False):
            try:
                with open(path, "r") as file:
                    lines = file.readlines()

                if not lines:
                    self.logger.error(f"File at {path} is empty.")
                    return []

                headers = lines[0].strip().split(",")
                if flag:
                    headers = headers[2:]

                data = []

                for line in lines[1:]:
                    self.logger.debug(f"Processing line: {line.strip()}")
                    values = line.strip().split(",")
                    values = values[2:] if flag else values

                    def infer_type(value): # TODO: Datastructure support
                        try:
                            return int(value)
                        except ValueError:
                            try:
                                return float(value)
                            except ValueError:
                                return value

                    _values = [infer_type(x) for x in values]

                    data.append(dict(zip(headers, _values)))

                return data
            except Exception as e:
                raise RuntimeError(
                    f"Failed to read CSV from {path}."
                ) from e


        def run_subprocess(command, env_vars=None):
            try:
                result = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    check=True,
                    env=env_vars,
                )
                self.logger.debug(
                    f"Subprocess output: {result.stdout}"
                )
                return result.stdout
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    "Subprocess execution failed."
                ) from e

        with tempfile.TemporaryDirectory() as temp_dir:

            model_path       = os.path.join(self.dir)
            temp_script_path = os.path.join(temp_dir, "script.sh")
            bash_output_path = os.path.join(temp_dir, "bash_output.csv")
            output_path      = os.path.join(temp_dir, "ersilia_output.csv")
            output_log_path  = os.path.join(temp_dir, "output.txt")
            error_log_path   = os.path.join(temp_dir, "error.txt")
            ex_file = os.path.join(temp_dir, "example_file.csv")

            self.run_exampe(
                n_samples=Options.NUM_SAMPLES.value,
                file_name=ex_file,
                simple=True,
                try_predefined=False
            )

            run_sh_path = os.path.join(
                model_path,
                "model",
                "framework",
                RUN_FILE
            )
            if not os.path.exists(run_sh_path):
                self.logger.warning(
                    f"{RUN_FILE} not found at {run_sh_path}. Skipping bash run."
                )
                return

            bash_script = f"""
                source {self._conda_prefix(self._is_base())}/etc/profile.d/conda.sh
                conda activate {self.model_id}
                cd {os.path.dirname(run_sh_path)}
                bash run.sh . {ex_file} {bash_output_path} > {output_log_path} 2> {error_log_path}
                conda deactivate
                """

            with open(temp_script_path, "w") as script_file:
                script_file.write(bash_script)

            self.logger.debug(f"Running bash script: {temp_script_path}")
            out = run_subprocess(["bash", temp_script_path])
            self.logger.info(out)
            bsh_data = read_csv(bash_output_path)
            self.logger.info(f"Bash Data:{bsh_data}")
            self.logger.debug("Serving the model after run.sh")
            run_subprocess(
                ["ersilia", "-v",
                 "serve", self.model_id,
                ]
            )
            self.logger.debug(
                "Running model for bash data consistency checking"
            )
            run_subprocess(
                ["ersilia", "-v",
                 "run",
                 "-i", ex_file,
                 "-o", output_path
                ]
            )
            data = read_csv(output_path, flag=True)

            status = compare_outputs(bsh_data, data)
            return status


    @staticmethod
    def _default_env():
        if "CONDA_DEFAULT_ENV" in os.environ:
            return os.environ["CONDA_DEFAULT_ENV"]
        return None

    @staticmethod
    def _conda_prefix(is_base):
        o = run_command_check_output("which conda").rstrip()
        if o:
            o = os.path.abspath(os.path.join(o, "..", ".."))
            return o
        if is_base:
            o = run_command_check_output("echo $CONDA_PREFIX").rstrip()
            return o
        else:
            o = run_command_check_output("echo $CONDA_PREFIX_1").rstrip()
            return o

    def _is_base(self):
        default_env = self._default_env()
        self.logger.debug(f"Default environment: {default_env}")
        return default_env == "base"

    def _compare_string_similarity(
            self,
            str1,
            str2,
            threshold
        ):
        similarity = fuzz.ratio(str1, str2)
        return similarity >= threshold

    def run(self):
        """
        Run the model tests and checks.

        Raises
        ------
        ImportError
            If required packages are missing.
        """
        try:
            echo("* Basic test started: downloading the model")
            if self.from_dockerhub:
                self.setup_service.from_github = True
            self.setup_service.get_model()
            self._perform_basic_checks()
            echo("* Basic checks done!")
            if self.shallow:
                echo("Performing shallow checks")
                self._perform_shallow_checks()
                echo("Shallow checks done")
            if self.deep:
                self._perform_shallow_checks()
                self._perform_deep_checks()

        except Exception as e:
            echo(
                f"An error occurred: {e}"
            )
        finally:
            echo(
                "Run process finished successfully."
            )

    def _perform_basic_checks(self):
        self.checkup_service.check_information()
        self._generate_table(
            **TABLE_CONFIGS[TableType.MODEL_INFORMATION_CHECKS].__dict__,
            rows=self.ios_service.check_results,
        )
        self.ios_service.check_results.clear()

        self.checkup_service.check_files()
        self._generate_table(
            **TABLE_CONFIGS[TableType.MODEL_FILE_CHECKS].__dict__,
            rows=self.ios_service.check_results
        )
        self._log_directory_sizes()
        self._docker_yml_check()
    
    @show_loader(text="Performing shallow checks", color="cyan")
    def _perform_shallow_checks(self):
        self.fetch()
        out = self.checkup_service.check_model_output_content(self.run_exampe, self.run_model)
        res = []
        if self.from_github or self.from_s3:
            out1 = self._log_env_sizes()
            res.append(*out1)
        if self.from_dockerhub:
            size = self.ios_service.calculate_image_size()
            res.append(("Docker Image Size", size))

        self.logger.debug("Running model with single input...")
        out2 = self.checkup_service.check_single_input(
            self.run_model,
            self.run_exampe
        )
        out3 = self.checkup_service.check_example_input(
            self.run_model,
            self.run_exampe
        )
        out4 = self.checkup_service.check_consistent_output(
            self.run_exampe,
            self.run_model
        )
        res.append(*out2)
        res.append(*out3)
        res.append(*out4)
        self.ios_service._generate_table(
            **TABLE_CONFIGS[TableType.SHALLOW_CHECK_SUMMARY].__dict__,
            rows=res,
            large_table=True,
        )
        out5 = self.run_bash()
        self.ios_service._generate_table(
            **TABLE_CONFIGS[TableType.CONSISTENCY_BASH].__dict__,
            rows=[*out5],
        )
        self.ios_service._generate_table(
            **TABLE_CONFIGS[TableType.MODEL_OUTPUT].__dict__,
            rows=out,
        )

    @show_loader(text="Performing deep checks", color="cyan")
    def _perform_deep_checks(self):
        data = self.inspecter.run(["computational_performance_tracking"])
        self.logger.info(f"Data: {data}")
        self.ios_service._generate_table(
            **TABLE_CONFIGS[TableType.COMPUTATIONAL_PERFORMANCE].__dict__,
            rows=data,
            large_table=True,
        )
    
    def _docker_yml_check(self):
        data = self.inspecter.run(["docker_check"])
        self.ios_service._generate_table(
            **TABLE_CONFIGS[TableType.DEPENDECY_CHECK].__dict__,
            rows=data,
            large_table=True,
        )

    def _log_env_sizes(self):
        env_size = self.ios_service.get_env_sizes()
        return [("Env Size", env_size)]

    def _log_directory_sizes(self):
        dir_size = self.ios_service.get_directories_sizes()
        self._generate_table(
            **TABLE_CONFIGS[TableType.MODEL_DIRECTORY_SIZES].__dict__,
            rows=[("Directory", dir_size)]
        )

    def _generate_table(self, title, headers, rows, large=False):
        self.ios_service._generate_table(
            title=title,
            headers=headers,
            rows=rows,
            large_table=large
        )

class ModelTester(ErsiliaBase):
    """
    Class to handle model testing. Initializes the model tester services and runs the tests.
    Parameters
    ----------
    model_id : str
        The ID of the model.
    level : str
        The level of testing.
    dir : str
        The directory for the model.
    inspect : bool
        Whether to inspect the model.
    remote : bool
        Whether to fetch the model from a remote source.
    remove : bool
        Whether to remove the model after testing.
    """
    def __init__(
            self,
            model, 
            level, 
            from_dir, 
            from_github, 
            from_dockerhub, 
            from_s3, 
            version, 
            shallow,
            deep,
        ):
        ErsiliaBase.__init__(
            self,
            config_json=None,
            credentials_json=None
        )
        self.model_id = model
        self.level = level
        self.from_dir = from_dir
        self.model_dir = os.path.join(EOS_TMP, self.model_id)
        self.dir = from_dir or self.model_dir
        self.from_github = from_github
        self.from_dockerhub = from_dockerhub
        self.from_s3 = from_s3
        self.version = version
        self.shallow = shallow
        self.deep = deep
        self._check_pedendency()
        self.setup_service = SetupService(
            self.model_id,
            self.dir,
            self.from_github,
            self.from_s3,
            self.logger,
        )
        self.ios = IOService(
            self.logger,
            self.model_id,
            self.dir
        )
        self.checks = CheckService(
            self.logger,
            self.model_id,
            self.dir,
            self.ios,
        )
        self.inspecter = InspectService(
            dir=self.dir,
            model=self.model_id,
            remote=True
        )
        self.runner = RunnerService(
            self.model_id,
            self.logger,
            self.ios,
            self.checks,
            self.setup_service,
            self.level,
            self.dir,
            self._model_path,
            self.from_github,
            self.from_s3,
            self.from_dockerhub,
            self.version,
            self.shallow,
            self.deep,
            self.inspecter
        )
    def _check_pedendency(self):
        if MISSING_PACKAGES:
            raise ImportError(
                "Missing packages required for testing. "
                "Please install test extras with 'pip install ersilia[test]'."
            )

    def run(self):
        """
        Run the model tester.
        """
        self.runner.run()
