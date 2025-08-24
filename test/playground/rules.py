import ast
import csv
import json
import os
import re
import requests
from enum import Enum
from pathlib import Path
from ersilia.utils.conda import SimpleConda
from ersilia.utils.docker import SimpleDocker
from ersilia.utils.logging import logger
from ersilia.cli import echo
from ersilia.default import (
    EOS,
    DOCKERHUB_ORG,
    SESSIONS_DIR,
    SESSION_JSON,
    DOCKER_INFO_FILE,
    API_SCHEMA_FILE,
    STATUS_JOSN,
    INFORMATION_FILE,
)
from ersilia.utils.session import get_session_dir
from ersilia.publish.test.services.checks import CheckService
from ersilia.publish.test.services.io import IOService
from ersilia.publish.test.services.constants import STATUS_CONFIGS


class ResponseName(Enum):
    DEST_FOLDER_EXIST = "Dest Folder Exists"
    DEST_FOLDER_HAS_CONTENT = "Dest folder has necessary files"
    DOCKER_IMAGE_EXIST = "Docker Image Exists"
    CONDA_ENV_EXISTS = "Conda venv exists"
    SESSION_REQUIRED_FILES_EXIST = "Session required files exist"
    SESSION_FILE_EXISTS = "Session file exists"
    SESSION_SERVICE_CLASS_CORRECT = "Session service class is correct"
    SESSION_URL_VALID = "Session url is valid"
    FILE_CONTENT_VALID = "file content is valid"
    FILE_CONTENT_NOT_VALID = "file content is not valid"
    REPO_FOLDER_EXIST = "Repo Folder Exists"
    REPO_FOLDER_NOT_EXIST = "Repo Folder Not Exists"
    DEST_FOLDER_NOT_EXIST = "Dest Folder Not Existss"
    CONTAINERS_REMOVED = "Containers Not Existss"
    DOCKER_IMAGE_NOT_EXIST = "Docker Image Not Existss"
    CONDA_ENV_NOT_EXISTS = "Conda venv not exists"
    SESSION_FILES_NOT_EXIST = "Session files not exist"
    TEST_COMMAND_CHECKS_FAILED = "Test command checks"
    CATALOG_JSON_CONTENT_VALID = "Catalog json content is valid"
    CATALOG_TXT_FILE_CONTENT_VALID = "Catalog txt file content is valid"
    CATALOG_LOCAL_MODEL_COUNT_VALID = (
        "Local catalog count and fetched  model count matching"
    )
    VALID_EXAMPLE_JSON_CONTENT = "Valid example json content"
    VALID_EXAMPLE_CSV_CONTENT = "Valid example csv content and content length"
    VALID_EXAMPLE_CSV_LENGTH = "Valid content length"
    VALID_CATALOG_CONTENT = "Valid catalog content"


def create_response(name, status, detail=""):
    return {
        "name": name.value if isinstance(name, Enum) else name,
        "status": status,
        "detail": detail,
    }


RULE_REGISTRY = {}


def get_configs(command):
    flags = command[1]
    model = flags[0]
    source = flags[1] if len(flags) > 1 else ""
    version = flags[2] if len(flags) > 2 else ""
    return flags, model, source, version


class CommandRule:
    def check(self, *args, **kwargs):
        raise NotImplementedError("Each rule must implement a check method.")


def register_rule(name):
    def decorator(cls):
        RULE_REGISTRY[name] = cls
        return cls

    return decorator


@register_rule("fetch_rule")
class FetchRule(CommandRule):
    required_files = {
        API_SCHEMA_FILE,
        STATUS_JOSN,
        INFORMATION_FILE,
    }

    sources = {
        "github": "GitHub",
        "s3": "Amazon S3",
        "dockerhub": "DockerHub",
        "local": "Local Repository",
    }

    def __init__(self):
        self.sc = SimpleConda()
        self.sd = SimpleDocker()
        self.dest = os.path.join(EOS, "dest")
        self.repos = os.path.join(EOS, "repository")
        self.checkups = []
        self._invalid_files = []
        self._missed_files = []

    def check(self, command, config, std_out):
        self.checkups.append(self._check_folders(command))
        img_chks = self._check_docker_image(command)
        if img_chks:
            self.checkups.append(img_chks)
        self.checkups.append(self._check_dest_file_content(command))
        env_chks = self._check_conda_env(command)
        if env_chks:
            self.checkups.append(env_chks)
        return self.checkups

    def _check_folders(self, command):
        model = get_configs(command)[1]
        path = os.path.join(self.dest, model)
        if not os.path.exists(path):
            return create_response(
                name=ResponseName.DEST_FOLDER_EXIST,
                status=False,
                detail=f"Folder does not exist at: {path}",
            )
        return create_response(name=ResponseName.DEST_FOLDER_EXIST, status=True)

    def _check_dest_file_content(self, command):
        model = get_configs(command)[1]
        source = get_configs(command)[2]

        if "dockerhub" in source:
            self.required_files.add(DOCKER_INFO_FILE)

        dest_folder = os.path.join(self.dest, model)
        missing_files = {
            file
            for file in self.required_files
            if not os.path.exists(os.path.join(dest_folder, file))
        }

        self._missed_files.extend(missing_files)

        if missing_files:
            return create_response(
                name=ResponseName.DEST_FOLDER_HAS_CONTENT,
                status=False,
                detail=f"Missing files in destination folder: {missing_files}",
            )

        self._check_api_schema_file(dest_folder)
        self._check_docker_info_file(dest_folder, source)
        self._check_status_file(dest_folder)

        if self._invalid_files:
            return create_response(
                name=ResponseName.DEST_FOLDER_HAS_CONTENT,
                status=False,
                detail=f"Invalid content in files: {self._invalid_files}",
            )

        return create_response(
            name=ResponseName.DEST_FOLDER_HAS_CONTENT,
            status=True,
        )

    def _check_api_schema_file(self, dest_folder):
        api_schema_file = os.path.join(dest_folder, API_SCHEMA_FILE)
        try:
            with open(api_schema_file, "r") as file:
                json.load(file)
        except Exception as e:
            self._invalid_files.append(API_SCHEMA_FILE)

    def _check_docker_info_file(self, dest_folder, source):
        docker_info_file = os.path.join(dest_folder, DOCKER_INFO_FILE)
        error_message = f"Docker status found to be false in: {DOCKER_INFO_FILE}"
        if not os.path.exists(docker_info_file):
            return

        try:
            with open(docker_info_file, "r") as file:
                content = json.load(file)
                status_true = content.get("docker_hub")
                if "dockerhub" in source and not status_true:
                    self._invalid_files.append(error_message)
                elif not status_true:
                    self._invalid_files.append(error_message)

        except Exception as e:
            self._invalid_files.append(error_message)

    def _check_status_file(self, dest_folder):
        status_file = os.path.join(dest_folder, STATUS_JOSN)
        try:
            with open(status_file, "r") as file:
                content = json.load(file)
                if not content.get("done", False):
                    self._invalid_files.append(STATUS_JOSN)
        except Exception as e:
            self._invalid_files.append(STATUS_JOSN)

    def _check_docker_image(self, command):
        _, model, source, version = get_configs(command)
        version = version if version else "latest"
        if "dockerhub" in source:
            exist = self.sd.exists(org=DOCKERHUB_ORG, img=model, tag=version)
            return create_response(name=ResponseName.DOCKER_IMAGE_EXIST, status=exist)

    def _check_conda_env(self, command):
        model = get_configs(command)[1]
        source = get_configs(command)[2]
        if "github" in source or "s3" in source:
            exist = self.sc.exists(model)
            return create_response(name=ResponseName.CONDA_ENV_EXISTS, status=exist)


@register_rule("serve_rule")
class ServeRule(CommandRule):
    required_files = {
        "console.log",
        "current.log",
        "_logs",
        "logs",
        SESSION_JSON,
    }

    def __init__(self):
        self.session_dir = SESSIONS_DIR
        self.checkups = []

    def check(self, command, config, std_out):
        self.checkups.append(self.extract_url(std_out))
        return self.checkups

    def extract_url(self, text):
        pattern = r"http://[a-zA-Z0-9.-]+:\d+"
        match = re.search(pattern, text)
        url = match.group(0)
        response = requests.get(url)
        status_ok = response.status_code == 200
        return create_response(name=ResponseName.SESSION_URL_VALID, status=status_ok)


@register_rule("run_rule")
class RunRule(CommandRule):
    def __init__(self):
        self.ios = IOService(logger=logger, model_id=None, dir=None)
        self.check_service = CheckService(
            logger=logger,
            model_id=None,
            dir=None,
            from_github=None,
            from_dockerhub=None,
            ios=self.ios,
        )

    def check(self, command, config, std_out):
        return self.check_contents(command, config)

    def check_contents(self, command, config):
        flag = get_configs(command)[0]
        output_files = flag[3]
        input_data = flag[1]
        inp_type = self._get_input_type(command, "-i")
        echo(
            f"Input type: {inp_type} | Input: {input_data} | Output: {output_files}",
            fg="yellow",
            bold=True,
        )
        self.check_service.original_input_list = (
            self.check_service._get_original_input_list(inp_type, input_data)
        )
        res = self.check_service.validate_file_content(output_files, inp_type)
        if res[-1] == str(STATUS_CONFIGS.PASSED):
            return [
                create_response(
                    name=f"{output_files} {ResponseName.FILE_CONTENT_VALID}",
                    status=True,
                )
            ]
        else:
            return [
                create_response(
                    name=f"{output_files} {ResponseName.FILE_CONTENT_NOT_VALID}",
                    status=False,
                    detail=res[1],
                )
            ]

    def _get_input_type(self, command, flag_key=None):
        info = command[-1]
        match = re.search(r"flag:\s*(\[[^\]]+\])", info)
        if not match:
            raise ValueError("No flag information found in the provided command tuple.")

        flag_list_str = match.group(1)
        try:
            flag_list = ast.literal_eval(flag_list_str)
        except Exception as e:
            raise ValueError(f"Error parsing flag list: {e}")

        if not isinstance(flag_list, list) or len(flag_list) % 2 != 0:
            raise ValueError(
                "Flag information is not a valid list of flag-value pairs."
            )

        flag_dict = {
            flag_list[i]: flag_list[i + 1] for i in range(0, len(flag_list), 2)
        }

        if flag_key is not None:
            return flag_dict.get(flag_key)
        return flag_dict


@register_rule("delete_rule")
class DeleteRule(CommandRule):
    def __init__(self):
        self.dest = os.path.join(EOS, "dest")
        self.repos = os.path.join(EOS, "repository")
        self.sd = SimpleDocker()
        self.sc = SimpleConda()
        self.checkups = []

    def check(self, command, config, std_out):
        self.checkups.append(self._check_containers_images_removed(command))
        self.checkups.append(self._check_dest_folder(command))
        self.checkups.append(self._check_repo_folder(command))
        self.checkups.append(self._check_docker_image(command))
        self.checkups.append(self._check_conda_env(command))
        return self.checkups

    def _check_repo_folder(self, command):
        model = get_configs(command)[1]
        path = os.path.join(self.repos, model)
        if os.path.exists(path):
            return create_response(
                name=ResponseName.REPO_FOLDER_NOT_EXIST,
                status=False,
                detail=f"Folder does not exist at: {path}",
            )
        return create_response(name=ResponseName.REPO_FOLDER_EXIST, status=True)

    def _check_dest_folder(self, command):
        model = get_configs(command)[1]
        path = os.path.join(self.dest, model)
        if os.path.exists(path):
            return create_response(
                name=ResponseName.DEST_FOLDER_NOT_EXIST,
                status=False,
                detail=f"Folder does not exist at: {path}",
            )
        return create_response(name=ResponseName.DEST_FOLDER_EXIST, status=True)

    def _check_containers_images_removed(self, command):
        model = get_configs(command)[1]
        source = get_configs(command)[2]
        version = get_configs(command)[-1]
        version = version if version else "latest"
        img_name = f"{DOCKERHUB_ORG}/{model}:{version}"
        if "dockerhub" in source:
            containers = self.sd.containers(only_run=False)
            if containers:
                exists = img_name in containers.values()
                return create_response(
                    name=ResponseName.CONTAINERS_REMOVED, status=not exists
                )
        return create_response(name=ResponseName.CONTAINERS_REMOVED, status=True)

    def _check_docker_image(self, command):
        _, model, source, version = get_configs(command)
        version = version if version else "latest"
        if "dockerhub" in source:
            exist = self.sd.exists(org=DOCKERHUB_ORG, img=model, tag=version)
            return create_response(
                name=ResponseName.DOCKER_IMAGE_NOT_EXIST, status=exist
            )
        return create_response(name=ResponseName.DOCKER_IMAGE_NOT_EXIST, status=True)

    def _check_conda_env(self, command):
        _, model, source, _ = get_configs(command)
        if self.sc.exists(model):
            return create_response(name=ResponseName.CONDA_ENV_NOT_EXISTS, status=False)
        return create_response(name=ResponseName.CONDA_ENV_NOT_EXISTS, status=True)


@register_rule("close_rule")
class CloseRule(CommandRule):
    def __init__(self):
        self.current_session_dir = get_session_dir()

    def check(self, command, config, std_out):
        model = config["models"]["single"]
        return self._check_files(model)

    def _check_files(self, model):
        files = [
            f"{model}.pid",
            SESSION_JSON,
        ]
        exists = all(
            os.path.exists(os.path.join(EOS, self.current_session_dir, file))
            for file in files
        )
        if exists:
            return [
                create_response(name=ResponseName.SESSION_FILES_NOT_EXIST, status=False)
            ]
        return [
            create_response(
                name=ResponseName.SESSION_FILES_NOT_EXIST,
                status=True,
                detail=f"Session files do not exist at: {self.current_session_dir}",
            )
        ]


@register_rule("catalog_rule")
class CatalogRule(CommandRule):
    def __init__(self):
        self.dest = os.path.join(EOS, "dest")
        self.checkups = []

    def check(self, command, config, std_out):
        validators = [
            lambda: self.validate_json_data(std_out, command),
            lambda: self.validate_catalog_counts(std_out, command),
            lambda: self.validate_txt_data(command),
        ]

        self.checkups.extend(filter(None, (validator() for validator in validators)))

        return self.checkups

    def validate_json_data(self, std_out, command):
        flags = get_configs(command)[0]
        if "--as-json" in flags and std_out != "":
            data = json.loads(std_out)
            for obj in data:
                for key, value in obj.items():
                    if value is None or value == "" or value == [] or value == {}:
                        return create_response(
                            name=ResponseName.CATALOG_JSON_CONTENT_VALID,
                            status=False,
                            detail=f"Validation failed for key '{key}' in object: {obj}",
                        )
            return create_response(
                name=ResponseName.CATALOG_JSON_CONTENT_VALID,
                status=True,
            )

    def validate_txt_data(self, command):
        flags = get_configs(command)[0]
        has_file = any(flag.endswith(".txt") for flag in flags)
        if has_file:
            file_path = flags[flags.index("--file_name") + 1]
            with open(file_path, "r", encoding="utf-8") as file:
                reader = csv.DictReader(file)
                for row in reader:
                    for key, value in row.items():
                        if value is None or value.strip() == "":
                            return create_response(
                                name=ResponseName.CATALOG_TXT_FILE_CONTENT_VALID,
                                status=False,
                                detail=f"Validation failed for key '{key}' in object: {row}",
                            )
            return create_response(
                name=ResponseName.CATALOG_JSON_CONTENT_VALID,
                status=True,
            )

    def validate_catalog_counts(self, std_out, command):
        flags = get_configs(command)[0]
        if "--hub" not in flags and "--as-json" in flags:
            data_count = len(json.loads(std_out))
            model_count = sum(
                os.path.isdir(os.path.join(self.dest, d)) for d in os.listdir(self.dest)
            )
            is_match = data_count == model_count
            if is_match:
                return create_response(
                    name=ResponseName.CATALOG_LOCAL_MODEL_COUNT_VALID,
                    status=True,
                )
            return create_response(
                name=ResponseName.CATALOG_LOCAL_MODEL_COUNT_VALID,
                status=False,
                detail="Number of models fetched and accessed by catalog command don't not match",
            )


@register_rule("example_rule")
class ExampleRule(CommandRule):
    def __init__(self):
        self.checkups = []

    def check(self, command, config, std_out):
        response = self.validate_json_data(std_out, command, config)
        return [response]

    def validate_json_data(self, std_out, command, config):
        flags = get_configs(command)[0]
        nsample = self.get_nsamples_flag_value(flags)
        has_predefined = self._get_predefined(flags)
        if (
            not any(
                flag.endswith(".csv") for flag in flags if not isinstance(flag, int)
            )
            and std_out != ""
        ):
            data = json.loads(std_out)
            for obj in data:
                for key, value in obj.items():
                    if value is None or value == "" or value == [] or value == {}:
                        return create_response(
                            name=ResponseName.VALID_EXAMPLE_JSON_CONTENT,
                            status=False,
                            detail=f"Validation failed for key '{key}' in object: {obj}",
                        )
            if len(data) != nsample and has_predefined:
                return create_response(
                    name=ResponseName.VALID_EXAMPLE_CSV_LENGTH,
                    status=False,
                    detail=f"Incorrect number of example unputs found in json. Expected: \
                        {nsample}, Found: {len(len(data))}",
                )
            return create_response(
                name=ResponseName.VALID_EXAMPLE_JSON_CONTENT,
                status=True,
            )
        else:
            example_file = self.get_file_flag_value(flags)
            if not Path(example_file).exists():
                return False, f"Example file not found: {example_file}"

            response = self._validate_csv_file(example_file, nsample, flags)
            return response

    def _validate_csv_file(self, example_file, nsample, flags):
        has_predefined = self._get_predefined(flags)
        try:
            with open(example_file, mode="r", newline="", encoding="utf-8") as csvfile:
                reader = csv.DictReader(csvfile)
                _row_len = 0
                for row_number, row in enumerate(reader, start=1):
                    _row_len += 1
                    for column, value in row.items():
                        if value is None or value.strip() == "":
                            return create_response(
                                name=ResponseName.VALID_EXAMPLE_JSON_CONTENT,
                                status=False,
                                detail=f"Row {row_number}, Column '{column}' contains null or empty data.",
                            )
            if _row_len != nsample and has_predefined:
                return create_response(
                    name=ResponseName.VALID_EXAMPLE_CSV_LENGTH,
                    status=False,
                    detail=f"Incorrect number of example inputs found in csv. \
                        Expected: {nsample}, Found: {_row_len}",
                )
            return create_response(
                name=ResponseName.VALID_EXAMPLE_CSV_CONTENT,
                status=True,
            )
        except Exception as e:
            return False, f"Error reading CSV file: {str(e)}"

    def get_file_flag_value(self, flags):
        try:
            file_flag_index = (
                flags.index("--file_name")
                if "--file_name" in flags
                else flags.index("-f")
                if "-f" in flags
                else None
            )
            return flags[file_flag_index + 1] if file_flag_index is not None else None
        except (IndexError, ValueError):
            return None

    def get_nsamples_flag_value(self, flags):
        try:
            nsamples_flag_index = (
                flags.index("--n_samples")
                if "--n_samples" in flags
                else flags.index("-n")
                if "-n" in flags
                else None
            )
            return (
                flags[nsamples_flag_index + 1]
                if nsamples_flag_index is not None
                else None
            )
        except (IndexError, ValueError):
            return None

    def _get_predefined(self, flag):
        return "-p" in flag or "--predefined" in flag


@register_rule("test_rule")
class TestRule(CommandRule):
    def __init__(self):
        self.current_session_dir = get_session_dir()

    def check(self, command, config, std_out):
        model = config["models"]["single"]
        return self._check_status(model)
        

    def _check_status(self, model):
        file = f"{model}-test.json"
        with open(file, "r") as f:
            data = json.load(f)
        failed_checks = []
        for _, item in data.items():
            for key, value in item.items():
                if key in "computational_performance_tracking_details":
                    continue
                if not value:
                    error = f"{' '.join(key.split('_')).capitalize()} is failed"
                    failed_checks.append(error)
                    echo(error)
        if failed_checks:
            return [
            create_response(
                name=ResponseName.TEST_COMMAND_CHECKS_FAILED,
                status=False,
                detail=failed_checks,
            )
        ]
        return [
            create_response(
                name=ResponseName.TEST_COMMAND_CHECKS_FAILED,
                status=True,
                detail=f"All test command checks passed",
            )
        ]

def get_rule(rule_name, *args, **kwargs):
    rule_class = RULE_REGISTRY.get(rule_name)
    if not rule_class:
        raise ValueError(f"Rule '{rule_name}' is not registered.")
    return rule_class().check(*args, **kwargs)
