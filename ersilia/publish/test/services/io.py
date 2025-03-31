import csv
import docker
import json
import os
import platform
import re
import warnings
import yaml
import warnings

warnings.filterwarnings("ignore", message="Using slow pure-python SequenceMatcher")

# ruff: noqa
MISSING_PACKAGES = False
try:
    from rich.console import Console
    from rich.table import Table
    from rich.text import Text
except ImportError:
    MISSING_PACKAGES = True
# ruff: enable
from typing import List
from .... import throw_ersilia_exception
from ....default import (
    METADATA_JSON_FILE,
    METADATA_YAML_FILE,
    PACK_METHOD_BENTOML,
    PACK_METHOD_FASTAPI,
    EOS_TMP,
    DOCKERHUB_ORG,
    PREDEFINED_EXAMPLE_OUTPUT_FILES,
    PREDEFINED_EXAMPLE_INPUT_FILES,
    INSTALL_YAML_FILE,
    DOCKERFILE_FILE,
)
from .constants import (
    Checks,
    STATUS_CONFIGS,
    ERSILIAPACK_BACK_FILES,
    ERSILIAPACK_FILES,
    BENTOML_FILES,
)
from .setup import SetupService
from ....hub.fetch.actions.template_resolver import TemplateResolver
from ....utils.conda import SimpleConda
from ....utils.exceptions_utils import test_exceptions as texc
from ....utils.terminal import run_command
from ....utils.logging import logger
from ....utils.docker import SimpleDocker, set_docker_host
from ....cli import echo


class IOService:
    """
    Service for handling input/output operations related to model testing.

    Parameters
    ----------
    logger : logging.Logger
        Logger for logging messages.
    model_id : str
        Identifier of the model.
    dir : str
        Directory where the model repository is located.

    Examples
    --------
    .. code-block:: python

        ios = IOService(
            logger=logger,
            model_id="model_id",
            dir="/path/to/dir",
        )
        ios.read_metadata()
    """

    def __init__(self, logger, model_id: str, dir: str):
        self.logger = logger
        self.model_id = model_id
        self.dir = dir
        self.model_size = 0
        self.console = Console()
        self.check_results = []
        self.simple_docker = SimpleDocker()
        self.conda = SimpleConda()

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
        resolver = TemplateResolver(model_id=model_id, repo_path=repo_path)
        if resolver.is_bentoml():
            return PACK_METHOD_BENTOML
        elif resolver.is_fastapi():
            return PACK_METHOD_FASTAPI
        else:
            return None

    @staticmethod
    def read_csv_values(file_path):
        """
        Reads a CSV file and returns its values excluding the header row.

        Args:
            file_path (str): The path to the CSV file.

        Returns:
            list: A list of rows (each row is a list of values) from the CSV file without the header.
        """
        try:
            with open(file_path, "r", newline="") as csvfile:
                reader = csv.reader(csvfile)
                next(reader, None)
                return [row for row in reader]
        except Exception as e:
            return []

    @staticmethod
    def _get_output_file_path(dir):
        for filename in PREDEFINED_EXAMPLE_OUTPUT_FILES:
            output_path = os.path.join(dir, filename)
            if os.path.exists(output_path):
                return output_path
        return None

    @staticmethod
    def _get_input_file_path(dir):
        for filename in PREDEFINED_EXAMPLE_INPUT_FILES:
            input_path = os.path.join(dir, filename)
            if os.path.exists(input_path):
                return input_path
        return None

    @staticmethod
    def _get_input_from_example_file(dir):
        input_path = IOService._get_input_file_path(dir)
        if input_path:
            inputs = IOService.read_csv_values(input_path)
            inputs = [input[0] for input in inputs]
            return inputs

    def get_output_consistency(self):
        data = self._read_metadata()
        if data is None:
            return "Fixed"
        if "Output Consistency" in data:
            return data["Output Consistency"]
        return "Fixed"

    def get_file_requirements(self):
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
        model_type = IOService.get_model_type(
            model_id=self.model_id, repo_path=self.dir
        )
        if model_type == PACK_METHOD_BENTOML:
            return BENTOML_FILES
        elif model_type == PACK_METHOD_FASTAPI:
            if os.path.exists(os.path.join(EOS_TMP, self.model_id, METADATA_JSON_FILE)):
                return ERSILIAPACK_BACK_FILES
            else:
                return ERSILIAPACK_FILES
        else:
            raise None

    def _run_check(
        self, check_function, data, check_name: str, additional_info=None
    ) -> bool:
        try:
            if additional_info is not None:
                check_function(additional_info)
            else:
                check_function(data)
            details = (
                f"Field {check_name} has correct entry." if data else "File exists"
            )
            self.check_results.append((check_name, details, str(STATUS_CONFIGS.PASSED)))
            return True

        except texc.EmptyField as ef:
            self.logger.error(f"EmptyField exception caught for key '{ef}'")
            details = (
                f"Field {check_name} has invalid value"
                if data
                else f"File {check_name} does not exist"
            )
            self.check_results.append((check_name, details, str(STATUS_CONFIGS.FAILED)))
            return False

        except texc.EmptyKey as ek:
            self.logger.error(f"EmptyKey exception caught for key '{ek}'")
            details = (
                f"Key {check_name} not present in the metadata file"
                if data
                else f"File {check_name} does not exist"
            )
            self.check_results.append(
                (check_name, details, str(STATUS_CONFIGS.NOT_PRESENT))
            )
            return False

        except Exception as e:
            self.logger.error(f"An unexpected exception occurred: {e}")
            details = (
                f"Field {check_name} has invalid value"
                if data
                else f"File {check_name} does not exist"
            )
            self.check_results.append((check_name, details, str(STATUS_CONFIGS.FAILED)))
            return False

    def _get_metadata_file(self):
        if METADATA_JSON_FILE in self.get_file_requirements():
            path = os.path.join(self.dir, METADATA_JSON_FILE)
        elif METADATA_YAML_FILE in self.get_file_requirements():
            path = os.path.join(self.dir, METADATA_YAML_FILE)
        else:
            return None
        return path

    def _read_metadata(self):
        path = self._get_metadata_file()
        if path is not None:
            with open(path, "r") as file:
                if path.endswith(".json"):
                    data = json.load(file)
                elif path.endswith((".yml", ".yaml")):
                    data = yaml.safe_load(file)
                else:
                    raise ValueError(f"Unsupported file format: {path}")
            return data
        return None

    def _combine_dir_size(self, data):
        keys = ("directory_size_mb", "model_size_check", "directory_size_check")
        if keys[1] in data:
            data[keys[1]][keys[0]] = data[keys[2]][keys[0]]
            del data[keys[2]]
            return data
        return data

    def collect_and_save_json(self, results, output_file):
        """
        Helper function to collect JSON results and save them to a file.
        """
        aggregated_json = {}
        for result in results:
            self.logger.info(f"Json result\n: {result}")
            aggregated_json.update(result)
        aggregated_json = self._combine_dir_size(aggregated_json)

        with open(output_file, "w") as f:
            json.dump(aggregated_json, f, indent=4)

    def _create_json_data(self, rows, key):
        def clean_string(s):
            s = s.replace("\n", "")
            s = re.sub(r"\\u[0-9a-fA-F]{4}", "", s)
            s = s.strip()
            return s

        def sanitize_name(name):
            return re.sub(r"[ \-./]", "_", str(name).lower())

        def parse_status(status):
            if isinstance(status, str):
                status = re.sub(
                    r"[-+]?\d*\.\d+|\d+", lambda m: str(float(m.group())), status
                )
                if "[green]✔" in status:
                    return True
                elif "[red]✘" in status:
                    return False
                elif "[yellow]✘" in status:
                    return "not present"
                else:
                    return re.sub(r"\[.*?\]", "", status).strip()
            return status

        def parse_performance(status):
            return {
                f"pred_{match[0]}": float(match[1])
                for match in re.findall(
                    r"(\d+) predictions executed in (\d+\.\d{2}) seconds. \n", status
                )
            }

        key = clean_string(sanitize_name(key))
        json_data = {}

        for row in rows:
            check_name = sanitize_name(row[0])
            check_status = row[-1]

            if check_name == "computational_performance_tracking_details":
                json_data[check_name] = parse_performance(check_status)
            elif check_name in (
                "environment_size_mb",
                "directory_size_mb",
                "image_size_mb",
            ):
                json_data[check_name] = float(check_status)
            else:
                json_data[check_name] = parse_status(clean_string(check_status))

        return {key: json_data}

    def _generate_table(self, title, headers, rows, large_table=False, merge=False):
        f_col_width = 30 if large_table else 30
        l_col_width = 50 if large_table else 10
        d_col_width = 30 if not large_table else 20

        table = Table(
            title=Text(title, style="bold light_green"),
            border_style="light_green",
            show_lines=True,
        )

        table.add_column(headers[0], justify="left", width=f_col_width, style="bold")
        for header in headers[1:-1]:
            table.add_column(header, justify="center", width=d_col_width, style="bold")
        table.add_column(headers[-1], justify="right", width=l_col_width, style="bold")

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
                row[-1],
            ]
            table.add_row(*styled_row)

        json_data = self._create_json_data(rows, title)

        self.console.print(table)

        return json_data

    def get_conda_env_size(self):
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
            loc = SetupService.get_conda_env_location(self.model_id, self.logger)
            return self.calculate_directory_size(loc)
        except Exception as e:
            self.logger.error(
                f"Error calculating size of Conda environment '{self.model_id}': {e}"
            )
            echo(f"Error calculating size of Conda environment '{self.model_id}': {e}")
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
                shell=False,
            )
            size = float(size_output.split()[0])
            return size
        except Exception as e:
            self.logger.error(f"Error calculating directory size for {path}: {e}")
            echo(f"Error calculating directory size for {path}: {e}")
            return 0

    def calculate_image_size(self, tag="latest"):
        """
        Calculate the size of a Docker image.

        Parameters
        ----------
        tag : str, optional
            The tag of the Docker image (default is 'latest').

        Returns
        -------
        str
            The size of the Docker image.
        """
        image_name = f"{DOCKERHUB_ORG}/{self.model_id}:{tag}"
        set_docker_host()
        client = docker.from_env()
        try:
            image = client.images.get(image_name)
            size_bytes = image.attrs["Size"]
            size_mb = size_bytes / (1024**2)
            return f"{size_mb:.2f}", Checks.SIZE_CACL_SUCCESS.value
        except docker.errors.ImageNotFound:
            return f"Image '{image_name}' not found.", Checks.SIZE_CACL_FAILED.value
        except Exception as e:
            return f"An error occurred: {e}", Checks.SIZE_CACL_FAILED.value

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
        dir_size = f"{dir_size:.2f}"
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
        env_size = f"{env_size:.2f}"
        return env_size


class PackageInstaller:
    def __init__(self, dir, model_id):
        self.dir = dir
        self.model_id = model_id
        self.conda = SimpleConda()
        self.logger = logger

    def _get_python_version_from_line(self, line: str):
        m = re.search(r"-py(\d+)", line)
        if m:
            ver_str = m.group(1)
            return f"{ver_str[0]}.{ver_str[1:]}"
        return None

    def _parse_install_yaml(self, path: str):
        with open(path, "r") as f:
            data = yaml.safe_load(f)
        python_version = data.get("python")
        pip_packages = []
        conda_packages = []
        commands = data.get("commands", [])
        for cmd in commands:
            if cmd[0].lower() == "pip":
                pkg_spec = f"{cmd[1]}=={cmd[2]}" if len(cmd) >= 3 else cmd[1]
                pip_packages.append(pkg_spec)
            elif cmd[0].lower() == "conda":
                pkg_spec = f"{cmd[1]}={cmd[2]}" if len(cmd) >= 3 else cmd[1]
                conda_packages.append(pkg_spec)
        return python_version, pip_packages, conda_packages

    def _parse_dockerfile(self, path: str):
        pip_packages = []
        conda_packages = []
        python_version = None
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("FROM"):
                    python_version = self._get_python_version_from_line(line)
                if line.startswith("RUN pip install"):
                    cmd = line.split("#")[0]
                    parts = cmd.split()
                    if len(parts) >= 4:
                        pkg_spec = parts[3]
                        pip_packages.append(pkg_spec)
                if line.startswith("RUN conda install"):
                    cmd = line.split("#")[0]
                    parts = cmd.split()
                    if parts:
                        pkg_spec = parts[-1]
                        conda_packages.append(pkg_spec)
        return python_version, pip_packages, conda_packages

    def _install_packages_from_dir(self):
        install_yml_path = os.path.join(self.dir, INSTALL_YAML_FILE)
        dockerfile_path = os.path.join(self.dir, DOCKERFILE_FILE)

        if os.path.exists(install_yml_path):
            python_version, pip_packages, conda_packages = self._parse_install_yaml(
                install_yml_path
            )
        elif os.path.exists(dockerfile_path):
            python_version, pip_packages, conda_packages = self._parse_dockerfile(
                dockerfile_path
            )
        else:
            self.logger.info(
                "Neither 'install.yml' nor 'Dockerfile' was found in the specified directory."
            )
            return
        
        if self.conda.exists(self.model_id):
            return
        
        if not self.conda.exists(self.model_id):
            self.conda.create(self.model_id, python_version)
        for pkg in pip_packages:
            self.logger.info(f"Installing pip package: {pkg}")
            run_command(["conda", "run", "-n", self.model_id, "pip", "install", pkg])
        for pkg in conda_packages:
            self.logger.info(f"Installing conda package: {pkg}")
            run_command(["conda", "install", "-n", self.model_id, "-y", pkg])
        self.logger.info(
            f"Installation complete in the conda environment: {self.model_id}"
        )
