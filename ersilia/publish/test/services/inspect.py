import os
import re
import subprocess
import time
from collections import namedtuple

import requests
import yaml

from ....default import (
    DOCKERFILE_FILE,
    INSTALL_YAML_FILE,
    METADATA_JSON_FILE,
    METADATA_YAML_FILE,
    PACK_METHOD_BENTOML,
    PACK_METHOD_FASTAPI,
    DOCKERFILE_FILE,
)
from .io import IOService
from .constants import STATUS_CONFIGS, Options
from .constants import (
    BASE_URL,
    RAW_CONTENT_URL,
    REPO_API_URL,
    USER_AGENT,
    COMMON_FILES,
    BENTOML_FILES,
    ERSILIAPACK_FILES,
)
from .... import ErsiliaBase
from ....hub.content.card import RepoMetadataFile
from ....hub.fetch.actions.template_resolver import TemplateResolver
from ....utils.logging import logger

Result = namedtuple("Result", ["success", "details"])


class ModelInspector:
    """
    Class for inspecting model repositories.

    Parameters
    ----------
    model : str
        The ID of the model to be inspected.
    dir : str
        The directory where the model repository is located.
    config_json : str, optional
        Path to the configuration JSON file.

    Examples
    --------
    .. code-block:: python

        inspector = ModelInspector(
            model="model_id", dir="path/to/repo"
        )
        result = inspector.check_repo_exists()
        result = inspector.check_complete_metadata()
    """

    BENTOML_FILES = COMMON_FILES + BENTOML_FILES
    ERSILIAPACK_FILES = COMMON_FILES + ERSILIAPACK_FILES

    REQUIRED_FIELDS = ["Publication", "Source Code", "S3", "DockerHub"]

    def __init__(self, model: str, dir: str, config_json=None):
        self.model = model
        self.dir = dir
        self.repo_url = f"{BASE_URL}{model}"
        self.content_url = RAW_CONTENT_URL.format(model=model)
        self.config_json = config_json
        self.pack_type = IOService.get_model_type(self.model, self.dir)

    def check_repo_exists(self):
        """
        Check if the model repository exists.

        Returns
        -------
        Result
            A namedtuple containing the success status and details of the check.
        """
        if self._url_exists(self.repo_url):
            return Result(True, "Repository exists.")
        return Result(False, f"Repository not found at {self.repo_url}.")

    def check_complete_metadata(self):
        """
        Check if the metadata file is complete.

        Returns
        -------
        Result
            A namedtuple containing the success status and details of the check.
        """
        url = (
            f"{self.content_url}{METADATA_JSON_FILE}"
            if self.pack_type == "bentoml"
            else f"{self.content_url}{METADATA_YAML_FILE}"
        )
        if not self._url_exists(url):
            return Result(False, f"Metadata file missing at {url}.")

        metadata = self._fetch_json(url)
        if metadata is None:
            return Result(False, "Failed to fetch or parse metadata.")

        missing_fields = [
            field for field in self.REQUIRED_FIELDS if field not in metadata
        ]
        invalid_urls = [
            (field, metadata[field])
            for field in self.REQUIRED_FIELDS
            if field in metadata and not self._url_exists(metadata[field])
        ]

        details = []
        if missing_fields:
            details.append(
                f"Missing fields: {', '.join(missing_fields)}. \
                  Required: {', '.join(self.REQUIRED_FIELDS)}."
            )
        if invalid_urls:
            details.extend(
                f"Invalid URL in '{field}': {url}" for field, url in invalid_urls
            )

        try:
            RepoMetadataFile.read_information(RepoMetadataFile(self.model))
        except Exception as e:
            details.append(f"Error encountered when parsing metadata file: {e}")

        if details:
            return Result(False, " ".join(details))

        return Result(True, "Metadata is complete.")

    def check_dependencies_are_valid(self):
        """
        Check if the dependencies in the Dockerfile or install.yml are valid.

        For PACK_METHOD_BENTOML, the Dockerfile is validated.
        For PACK_METHOD_FASTAPI, either file is acceptable. If one exists, it is validated;
        if both exist, both are validated. If neither is present, an error is returned.

        Returns
        -------
        Result
            A namedtuple containing the success status and details of the check.
        """
        if self.pack_type not in [PACK_METHOD_BENTOML, PACK_METHOD_FASTAPI]:
            return Result(False, f"Unsupported pack type: {self.pack_type}")

        if self.pack_type == PACK_METHOD_BENTOML:
            content, error = self._get_file_content(DOCKERFILE_FILE)
            if content is None:
                return Result(False, error)
            errors = self._validate_dockerfile(content)
            if errors:
                return Result(False, " ".join(errors))
            return Result(True, f"{DOCKERFILE_FILE} dependencies are valid.")

        else:
            dockerfile_content, dockerfile_error = self._get_file_content(
                DOCKERFILE_FILE
            )
            yml_content, yml_error = self._get_file_content(INSTALL_YAML_FILE)

            if dockerfile_content is None and yml_content is None:
                return Result(
                    False,
                    f"Neither {DOCKERFILE_FILE} nor {INSTALL_YAML_FILE} found. "
                    f"Dockerfile error: {dockerfile_error} | "
                    f"install.yml error: {yml_error}",
                )

            errors = []

            if dockerfile_content is not None:
                dockerfile_errors = self._validate_dockerfile(dockerfile_content)
                if dockerfile_errors:
                    errors.extend(dockerfile_errors)

            if yml_content is not None:
                yml_errors = self._validate_yml(yml_content)
                if yml_errors:
                    errors.extend(yml_errors)

            if errors:
                return Result(False, " ".join(errors))

            valid_files = []
            if dockerfile_content is not None:
                valid_files.append(DOCKERFILE_FILE)
            if yml_content is not None:
                valid_files.append(INSTALL_YAML_FILE)

            return Result(True, f"{', '.join(valid_files)} dependencies are valid.")

    def _get_file_content(self, file):
        if self.dir is not None:
            path = os.path.join(self.dir, file)
            if not os.path.isfile(path):
                return None, f"{file} not found at {path}"
            try:
                with open(path, "r") as file:
                    return file.read(), None
            except Exception as e:
                return (None, f"Failed to read {file} content: {str(e)}")
        else:
            url = f"{self.content_url}{file}"
            if not self._url_exists(url):
                return None, f"{file} not found at {url}"
            content = self._fetch_text(url)
            if content is None:
                return (None, f"Failed to fetch {file} content.")
            return content, None

    def check_complete_folder_structure(self):
        """
        Check if the folder structure of the repository is complete.

        Returns
        -------
        Result
            A namedtuple containing the success status and details of the check.
        """
        invalid_items = self.validate_repo_structure()
        if invalid_items:
            return Result(False, f"Missing folders: {', '.join(invalid_items)}")
        return Result(True, "Folder structure is complete.")

    def check_computational_performance(self):
        """
        Check the computational performance of the model.

        Returns
        -------
        Result
            A namedtuple containing the success status and details of the check.
        """
        details = []
        for n in (1, 10, 100):
            result = self._run_performance_check(n)
            if not result.success:
                return result
            details.append(result.details)
        return Result(True, " ".join(details))

    def check_no_extra_files(self):
        """
        Check if there are no extra files in the repository.

        Returns
        -------
        Result
            A namedtuple containing the success status and details of the check.
        """
        if self.pack_type == PACK_METHOD_BENTOML:
            expected_items = self.BENTOML_FILES + self.BENTOML_FOLDERS
        elif self.pack_type == PACK_METHOD_FASTAPI:
            expected_items = self.ERSILIAPACK_FILES + self.ERSILIAPACK_FOLDERS
        else:
            return Result(False, f"Unsupported pack type: {self.pack_type}")

        if self.dir is not None:
            unexpected_items = []
            for root, dirs, files in os.walk(self.dir):
                relative_path = os.path.relpath(root, self.dir)
                items_in_dir = [
                    os.path.join(relative_path, item) for item in files + dirs
                ]
                unexpected_items.extend(
                    item for item in items_in_dir if item not in expected_items
                )

            if unexpected_items:
                return Result(
                    False, f"Unexpected items found: {', '.join(unexpected_items)}"
                )
            return Result(True, "No extra files found locally.")

        url = REPO_API_URL.format(model=self.model)
        if not self._url_exists(url):
            return Result(False, f"Failed to access repository contents at: {url}")
        headers = {
            "Accept": "application/vnd.github.v3+json",
        }

        response = requests.get(url, headers=headers)
        if response.status_code != 200:
            return Result(False, "Failed to fetch repository contents.")
        unexpected_items = [
            item["name"]
            for item in response.json()
            if item["name"] not in expected_items
        ]

        if unexpected_items:
            return Result(
                False, f"Unexpected items found: {', '.join(unexpected_items)}"
            )

        return Result(True, "No extra files found.")

    def _url_exists(self, url):
        try:
            headers = {
                "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
                "User-Agent": USER_AGENT,
            }
            response = requests.head(url, headers=headers)
            logger.debug(f"URl: {url} | status code: {response.status_code}")
            return response.status_code == 200
        except requests.RequestException:
            return False

    def _fetch_json(self, url):
        try:
            response = requests.get(url)
            return response.json()
        except (requests.RequestException, ValueError):
            return None

    def _fetch_text(self, url):
        try:
            response = requests.get(url)
            return response.text
        except requests.RequestException:
            return None

    def _validate_urls(self, metadata, fields):
        invalid_urls = []
        for field in fields:
            url = metadata.get(field)
            if url and not self._url_exists(url):
                invalid_urls.append((field, url, 404))
        return invalid_urls

    def _validate_repo_structure(self, required_items):
        missing_items = []

        if self.dir is not None:
            for item in required_items:
                item_path = os.path.join(self.dir, item)
                if not os.path.isfile(item_path):
                    missing_items.append(item)
        else:
            for item in required_items:
                url = f"{RAW_CONTENT_URL.format(model=self.model)}{item}"
                response = requests.head(url)
                if response.status_code != 200:
                    logger.debug(f"URL: {url} | STatus Code: {response.status_code}")
                    missing_items.append(item)

        return missing_items

    def validate_repo_structure(self):
        """
        Validate the repository structure.

        Returns
        -------
        List of missing items.
        """
        logger.debug(f"Pack Type: {self.pack_type}")
        if self.pack_type == PACK_METHOD_BENTOML:
            required_items = self.BENTOML_FILES
        elif self.pack_type == PACK_METHOD_FASTAPI:
            required_items = self.ERSILIAPACK_FILES
        else:
            raise ValueError(f"Unsupported pack type: {self.pack_type}")

        return self._validate_repo_structure(required_items)

    def _validate_dockerfile(self, dockerfile_content):
        lines, errors = dockerfile_content.splitlines(), []

        if "WORKDIR /repo" not in dockerfile_content:
            errors.append("Missing 'WORKDIR /repo'.")
        if "COPY . /repo" not in dockerfile_content:
            errors.append("Missing 'COPY . /repo'.")

        pip_install_pattern = re.compile(r"pip install (.+)")
        version_pin_pattern = re.compile(
            r"^[a-zA-Z0-9_\-\.]+(==|>=|<=|>|<)[a-zA-Z0-9_\-\.]+$"
        )

        for line in lines:
            line = line.strip()

            match = pip_install_pattern.search(line)
            if match:
                install_cmd, _, _ = match.group(1).partition(
                    "#"
                )  # Remove comments, this actually crashes things
                packages_and_options = install_cmd.strip().split()
                skip_next = False

                for item in packages_and_options:
                    if skip_next:
                        skip_next = False
                        continue

                    if item.startswith("--index-url") or item.startswith(
                        "--extra-index-url"
                    ):
                        skip_next = True
                        continue

                    if item.startswith("git+"):
                        continue

                    if not version_pin_pattern.match(item):
                        errors.append(
                            f"Package '{item}' in line '{line}' is not version-pinned (e.g., 'package==1.0.0')."
                        )

        return errors

    def _validate_yml(self, yml_content):
        errors = []
        try:
            yml_data = yaml.safe_load(yml_content)
        except yaml.YAMLError as e:
            return [f"YAML parsing error: {str(e)}"]

        python_version = yml_data.get("python")
        if not python_version:
            errors.append("Missing Python version in dependency.yml.")
        elif not re.match(r"^\d+(\.\d+){0,2}$", python_version):
            errors.append(f"Invalid Python version format: {python_version}")

        version_pin_pattern = re.compile(
            r"^[a-zA-Z0-9_\-\.]+(==|>=|<=|>|<)[a-zA-Z0-9_\-\.]+$"
        )

        commands = yml_data.get("commands", [])
        for command in commands:
            if not isinstance(command, list) or len(command) < 2:
                errors.append(f"Invalid command format: {command}")
                continue

            tool = command[0]
            package = command[1]
            version = command[2] if len(command) > 2 else None
            additional_args = command[3:] if len(command) > 3 else []
            if tool not in ("pip", "conda"):
                errors.append(f"Unsupported package manager: {tool}")
                continue
            if not version:
                errors.append(
                    f"Missing version for package '{package}' in command '{command}'."
                )
            if tool in ("pip", "conda"):
                if tool == "pip":
                    if package.startswith("git+"):
                        continue
                    if version and not version_pin_pattern.match(
                        f"{package}=={version}"
                    ):
                        errors.append(
                            f"Package '{package}' in command '{command}' is not properly version-pinned."
                        )

                    for arg in additional_args:
                        if arg.startswith("--index-url") or arg.startswith(
                            "--extra-index-url"
                        ):
                            continue
                        if not version_pin_pattern.match(arg):
                            errors.append(
                                f"Unexpected argument in command '{command}': {arg}"
                            )

                elif tool == "conda" and not version:
                    errors.append(
                        f"Package '{package}' in command '{command}' does not have a valid pinned version."
                    )

        return errors

    def _run_performance_check(self, n):
        cmd = (
            f"ersilia serve {self.model}&& "
            f"ersilia example -n {n} -c -f {Options.DEEP_INPUT.value} && "
            f"ersilia run -i {Options.DEEP_INPUT.value} && ersilia close"
        )
        start_time = time.time()
        process = subprocess.run(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if process.returncode != 0:
            return Result(False, f"Error serving model: {process.stderr.strip()}")
        execution_time = time.time() - start_time
        return Result(
            True, f"{n} predictions executed in {execution_time:.2f} seconds. \n"
        )


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
    remote : bool, optional
        Flag indicating whether the model is remote.
    config_json : str, optional
        Path to the configuration JSON file.
    credentials_json : str, optional
        Path to the credentials JSON file.

    Examples
    --------
    .. code-block:: python

        inspector = InspectService(
            dir="/path/to/model", model="model_id"
        )
        results = inspector.run()
    """

    def __init__(
        self,
        dir: str,
        model: str,
        remote: bool = False,
        config_json: str = None,
        credentials_json: str = None,
    ):
        super().__init__(config_json, credentials_json)
        self.dir = dir

        self.model = model
        self.remote = remote
        self.resolver = TemplateResolver(model_id=model, repo_path=self.dir)

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
            " ".join(word.capitalize() for word in k.split("_")): _transform_key(v)
            for k, v in output.items()
        }

        output = [(key, value) for key, value in output.items()]
        return output

    def _get_checks(self, inspector: ModelInspector) -> dict:
        def create_check(check_fn, key, details):
            return lambda: CheckStrategy(check_fn, key, details)

        docker_file_exists = os.path.isfile(os.path.join(self.dir, DOCKERFILE_FILE))
        dependency_check = "Dockerfile" if docker_file_exists else "Install_YAML"
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
