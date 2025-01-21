import os
import re
import subprocess
import time
from collections import namedtuple

import requests
import yaml

from ..default import (
    DOCKERFILE_FILE,
    INSTALL_YAML_FILE,
    METADATA_JSON_FILE,
    METADATA_YAML_FILE,
    PACK_METHOD_BENTOML,
    PACK_METHOD_FASTAPI,
    PREDEFINED_EXAMPLE_FILES,
    RUN_FILE,
)
from ..hub.content.card import RepoMetadataFile
from ..hub.fetch.actions.template_resolver import TemplateResolver
from ..utils.logging import logger

Result = namedtuple("Result", ["success", "details"])

# Base URL for the Ersilia OS Github
BASE_URL = "https://github.com/ersilia-os/"
RAW_CONTENT_URL = "https://raw.githubusercontent.com/ersilia-os/{model}/main/"
REPO_API_URL = "https://api.github.com/repos/ersilia-os/{model}/contents"
USER_AGENT = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36"


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

    RUN_FILE = f"model/framework/{RUN_FILE}"

    COMMON_FILES = [
        RUN_FILE,
        "README.md",
        "LICENSE",
    ]
    BENTOML_FOLDERS = ["model", "src", ".github"]
    BENTOML_FILES = [
        DOCKERFILE_FILE,
        METADATA_JSON_FILE,
        "src/service.py",
        "pack.py",
        ".gitignore",
        "input.csv",
    ]

    ERSILIAPACK_FOLDERS = ["model", ".github"]

    ERSILIAPACK_FILES = [
        INSTALL_YAML_FILE,
        METADATA_YAML_FILE,
        PREDEFINED_EXAMPLE_FILES[0],
        PREDEFINED_EXAMPLE_FILES[1],
        ".dockerignore",
        ".gitignore",
        ".gitattributes",
    ]

    BENTOML_FILES = COMMON_FILES + BENTOML_FILES
    ERSILIAPACK_FILES = COMMON_FILES + ERSILIAPACK_FILES

    REQUIRED_FIELDS = ["Publication", "Source Code", "S3", "DockerHub"]

    def __init__(self, model: str, dir: str, config_json=None):
        self.model = model
        self.dir = dir
        self.repo_url = f"{BASE_URL}{model}"
        self.content_url = RAW_CONTENT_URL.format(model=model)
        self.config_json = config_json
        self.pack_type = self.get_pack_type()

    def get_pack_type(self):
        """
        Determine the packaging method of the model.

        Returns
        -------
        str
            The packaging method, either 'bentoml' or 'fastapi'.
        """
        resolver = TemplateResolver(model_id=self.model, repo_path=self.dir)
        if resolver.is_bentoml():
            return PACK_METHOD_BENTOML
        elif resolver.is_fastapi():
            return PACK_METHOD_FASTAPI
        else:
            return None

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

        Returns
        -------
        Result
            A namedtuple containing the success status and details of the check.
        """
        if self.pack_type not in [PACK_METHOD_BENTOML, PACK_METHOD_FASTAPI]:
            return Result(False, f"Unsupported pack type: {self.pack_type}")

        file = (
            DOCKERFILE_FILE
            if self.pack_type == PACK_METHOD_BENTOML
            else INSTALL_YAML_FILE
        )
        method = (
            self._validate_dockerfile
            if self.pack_type == PACK_METHOD_BENTOML
            else self._validate_yml
        )

        content, error = self._get_file_content(file)
        if content is None:
            return Result(False, error)

        errors = method(content)
        if errors:
            return Result(False, " ".join(errors))

        return Result(True, f"{file} dependencies are valid.")

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
        version_pin_pattern = re.compile(r"^[a-zA-Z0-9_\-\.]+==[a-zA-Z0-9_\-\.]+$")

        for line in lines:
            line = line.strip()

            match = pip_install_pattern.search(line)
            if match:
                packages_and_options = match.group(1).split()
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
            errors.append("Missing Python version in install.yml.")

        version_pin_pattern = re.compile(r"^[a-zA-Z0-9_\-\.]+==[a-zA-Z0-9_\-\.]+$")

        commands = yml_data.get("commands", [])
        for command in commands:
            if not isinstance(command, list) or len(command) < 2:
                errors.append(f"Invalid command format: {command}")
                continue

            tool = command[0]
            _ = command[1]
            version = command[2] if len(command) > 2 else None

            if tool in ("pip", "conda"):
                if tool == "pip":
                    pip_args = command[1:]
                    skip_next = False

                    for item in pip_args:
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
                                f"Package '{item}' in command '{command}' is not version-pinned (e.g., 'package==1.0.0')."
                            )

                elif tool == "conda" and not version:
                    errors.append(
                        f"Package in command '{command}' does not have a valid pinned version "
                        f"(should be in the format ['conda', 'package_name', 'x.y.z'])."
                    )

        return errors

    def _run_performance_check(self, n):
        cmd = (
            f"ersilia serve {self.model}&& "
            f"ersilia example -n {n} -c -f my_input.csv && "
            "ersilia run -i my_input.csv && ersilia close"
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
