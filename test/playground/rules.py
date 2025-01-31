import csv
import json
import os
import re
import requests
import yaml
from enum import Enum
from pathlib import Path
from ersilia.utils.conda import SimpleConda
from ersilia.utils.docker import SimpleDocker
from ersilia.utils.logging import logger
from ersilia.default import (
    EOS,
    DOCKERHUB_ORG,
    SESSIONS_DIR,
    SESSION_JSON,
    MODEL_SOURCE_FILE,
    DOCKER_INFO_FILE,
    API_SCHEMA_FILE,
    STATUS_JOSN,
    INFORMATION_FILE,
)
from ersilia.utils.session import get_session_dir
from ersilia.publish.test import CheckService, IOService, STATUS_CONFIGS


class ResponseName(Enum):
    DEST_FOLDER_EXIST = "Dest Folder Exist"
    DEST_FOLDER_HAS_CONTENT = "Dest folder has necessary files"
    DOCKER_IMAGE_EXIST = "Docker Image Exist"
    CONDA_ENV_EXISTS = "Conda venv exists"
    SESSION_REQUIRED_FILES_EXIST = "Session required files exist"
    SESSION_FILE_EXISTS = "Session file exists"
    SESSION_SERVICE_CLASS_CORRECT = "Session service class is correct"
    SESSION_URL_VALID = "Session url is valid"
    FILE_CONTENT_VALID = "file content is valid"
    FILE_CONTENT_NOT_VALID = "file content is not valid"
    REPO_FOLDER_EXIST = "Repo Folder Exist"
    REPO_FOLDER_NOT_EXIST = "Repo Folder Not Exist"
    DEST_FOLDER_NOT_EXIST = "Dest Folder Not Exist"
    CONTAINERS_REMOVED = "Containers Removed"
    DOCKER_IMAGE_NOT_EXIST = "Docker Image Not Exist"
    CONDA_ENV_NOT_EXISTS = "Conda venv not exists"
    SESSION_FILES_NOT_EXIST = "Session files not exist"
    CATALOG_JSON_CONTENT_VALID = "Catalog json content is valid"
    VALID_EXAMPLE_JSON_CONTENT = "Valid example json content"
    VALID_EXAMPLE_CSV_CONTENT = (
        "Valid example csv content and content length"
    )
    VALID_EXAMPLE_CSV_LENGTH = "Valid content length"
    VALID_CATALOG_CONTENT = "Valid catalog content"


def create_response(name, status, detail=""):
    return {
        "name": name.value if isinstance(name, Enum) else name,
        "status": status,
        "detail": detail,
    }


RULE_REGISTRY = {}

config_path = Path("config.yml")
config = yaml.safe_load(config_path.read_text())


def get_configs(command):
    flags = command[1]
    model = flags[0]
    source = flags[1] if len(flags) > 1 else ""
    version = flags[2] if len(flags) > 2 else ""
    return flags, model, source, version


class CommandRule:
    def check(self, *args, **kwargs):
        raise NotImplementedError(
            "Each rule must implement a check method."
        )


def register_rule(name):
    def decorator(cls):
        RULE_REGISTRY[name] = cls
        return cls

    return decorator


@register_rule("fetch_rule")
class FetchRule(CommandRule):
    required_files = {
        MODEL_SOURCE_FILE,
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
        return create_response(
            name=ResponseName.DEST_FOLDER_EXIST, status=True
        )

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

        self._check_model_source_file(dest_folder, source)
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

    def _check_model_source_file(self, dest_folder, source):
        source = source.replace("--from_", "")
        if source not in self.sources:
            self._invalid_files.append(MODEL_SOURCE_FILE)
            return

        expected_value = self.sources[source]

        model_source_file = os.path.join(dest_folder, MODEL_SOURCE_FILE)
        try:
            with open(model_source_file, "r") as file:
                content = file.read().strip()
                if content != expected_value:
                    self._invalid_files.append(MODEL_SOURCE_FILE)
        except Exception as e:
            self._invalid_files.append(MODEL_SOURCE_FILE)

    def _check_api_schema_file(self, dest_folder):
        api_schema_file = os.path.join(dest_folder, API_SCHEMA_FILE)
        try:
            with open(api_schema_file, "r") as file:
                json.load(file)
        except Exception as e:
            self._invalid_files.append(API_SCHEMA_FILE)

    def _check_docker_info_file(self, dest_folder, source):
        docker_info_file = os.path.join(dest_folder, DOCKER_INFO_FILE)
        if not os.path.exists(docker_info_file):
            return

        try:
            with open(docker_info_file, "r") as file:
                content = json.load(file)
                if "dockerhub" in source and not content.get(
                    "docker_hub", False
                ):
                    self._invalid_files.append(DOCKER_INFO_FILE)
                elif "dockerhub" not in source and content.get(
                    "docker_hub", False
                ):
                    self._invalid_files.append(DOCKER_INFO_FILE)
        except Exception as e:
            self._invalid_files.append(DOCKER_INFO_FILE)

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
            exist = self.sd.exists(
                org=DOCKERHUB_ORG, img=model, tag=version
            )
            return create_response(
                name=ResponseName.DOCKER_IMAGE_EXIST, 
                status=exist
            )

    def _check_conda_env(self, command):
        model = get_configs(command)[1]
        source = get_configs(command)[2]
        if "github" in source or "s3" in source:
            exist = self.sc.exists(model)
            return create_response(
                name=ResponseName.CONDA_ENV_EXISTS, 
                status=exist
            )


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
        self.current_session_dir = get_session_dir()

    def check(self, command, config, std_out):
        runner = config["runtime"]["runner"]
        resp_one = self.get_latest_folder_and_check_files(command, runner)
        if resp_one:
            self.checkups.append(resp_one)
        resp_two = self.check_service_class(runner)
        if resp_two:
            self.checkups.append(resp_two)
        self.checkups.append(self.extract_url(std_out))
        return self.checkups

    def get_latest_folder_and_check_files(self, command, runner):
        if runner != "multiple":
            model = get_configs(command)[1]
            self.required_files.add(f"{model}.pid")
            existing_files = set(os.listdir(self.current_session_dir))
            required_files_exists = self.required_files.issubset(
                existing_files
            )
            if required_files_exists:
                return create_response(
                    name=ResponseName.SESSION_REQUIRED_FILES_EXIST,
                    status=True,
                )
            return create_response(
                name=ResponseName.SESSION_REQUIRED_FILES_EXIST,
                status=False,
                detail=f"Session required files do not exist at: {self.current_session_dir}",
            )

    def check_service_class(self, runner):
        if runner != "multiple":
            service_class = ("pulled_docker", "conda")
            session_file_path = os.path.join(
                self.current_session_dir, SESSION_JSON
            )

            if not os.path.exists(session_file_path):
                return create_response(
                    name=ResponseName.SESSION_FILE_EXISTS,
                    status=False,
                    detail=f"Session file does not exist at: {session_file_path}",
                )

            try:
                with open(session_file_path, "r") as file:
                    data = json.load(file)
                    correct_srv_class = (
                        data.get("service_class") in service_class
                    )
                    return create_response(
                        name=ResponseName.SESSION_SERVICE_CLASS_CORRECT,
                        status=correct_srv_class,
                        detail=f"Service class: {data.get('service_class')}",
                    )
            except (json.JSONDecodeError, KeyError):
                return create_response(
                    name=ResponseName.SESSION_SERVICE_CLASS_CORRECT,
                    status=False,
                    detail=f"Session file is not valid at: {session_file_path}",
                )

    def extract_url(self, text):
        pattern = r"http://[a-zA-Z0-9.-]+:\d+"
        match = re.search(pattern, text)
        url = match.group(0)
        response = requests.get(url)
        status_ok = response.status_code == 200
        return create_response(
            name=ResponseName.SESSION_URL_VALID, status=status_ok
        )


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
        inp_type = self._get_input_type(output_files)
        res = self.check_service.validate_file_content(
            output_files, inp_type
        )
        if res[-1] == str(STATUS_CONFIGS.PASSED):
            return [
                create_response(
                    name=ResponseName.FILE_CONTENT_VALID, status=True
                )
            ]
        else:
            return [
                create_response(
                    name=ResponseName.FILE_CONTENT_NOT_VALID,
                    status=False,
                    detail=res[1],
                )
            ]

    def _get_input_type(self, inp):
        if isinstance(inp, str):
            return "str"
        elif isinstance(inp, list):
            return "list"
        elif inp.endswith(".csv"):
            return "csv"
        return None


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
        print(self.checkups)
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
        return create_response(
            name=ResponseName.REPO_FOLDER_EXIST, status=True
        )

    def _check_dest_folder(self, command):
        model = get_configs(command)[1]
        path = os.path.join(self.dest, model)
        if os.path.exists(path):
            return create_response(
                name=ResponseName.DEST_FOLDER_NOT_EXIST,
                status=False,
                detail=f"Folder does not exist at: {path}",
            )
        return create_response(
            name=ResponseName.DEST_FOLDER_EXIST, status=True
        )

    def _check_containers_images_removed(self, command):
        model = get_configs(command)[1]
        version = get_configs(command)[-1]
        version = version if version else "latest"
        img_name = f"{DOCKERHUB_ORG}/{model}:{version}"
        containers = self.sd.containers(only_run=False)
        if containers:
            exists = img_name in containers.values()
            return create_response(
                name=ResponseName.CONTAINERS_REMOVED, status=not exists
            )
        return create_response(
            name=ResponseName.CONTAINERS_REMOVED, status=True
        )

    def _check_docker_image(self, command):
        _, model, source, version = get_configs(command)
        version = version if version else "latest"
        if "dockerhub" in source:
            exist = self.sd.exists(
                org=DOCKERHUB_ORG, img=model, tag=version
            )
            return create_response(
                name=ResponseName.DOCKER_IMAGE_NOT_EXIST, status=exist
            )
        return create_response(
            name=ResponseName.DOCKER_IMAGE_NOT_EXIST, status=True
        )

    def _check_conda_env(self, command):
        _, model, source, _ = get_configs(command)
        if self.sc.exists(model):
            return create_response(
                name=ResponseName.CONDA_ENV_NOT_EXISTS, status=False
            )
        return create_response(
            name=ResponseName.CONDA_ENV_NOT_EXISTS, status=True
        )


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
            os.path.exists(
                os.path.join(EOS, self.current_session_dir, file)
            )
            for file in files
        )
        if exists:
            return [
                create_response(
                    name=ResponseName.SESSION_FILES_NOT_EXIST, status=False
                )
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
        pass

    def check(self, command, config, std_out):
        response = self.validate_json_data(std_out, command)
        return [response]

    def validate_json_data(self, std_out, command):
        flags = get_configs(command)[0]
        if "--as-json" in flags:
            data = json.loads(std_out)
            for obj in data:
                for key, value in obj.items():
                    if (
                        value is None
                        or value == ""
                        or value == []
                        or value == {}
                    ):
                        return create_response(
                            name=ResponseName.CATALOG_JSON_CONTENT_VALID,
                            status=False,
                            detail=f"Validation failed for key '{key}' in object: {obj}",
                        )
            return create_response(
                name=ResponseName.CATALOG_JSON_CONTENT_VALID,
                status=True,
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

        if not any(
            flag.endswith(".csv")
            for flag in flags
            if not isinstance(flag, int)
        ):
            data = json.loads(std_out)
            for obj in data:
                for key, value in obj.items():
                    if (
                        value is None
                        or value == ""
                        or value == []
                        or value == {}
                    ):
                        return create_response(
                            name=ResponseName.VALID_EXAMPLE_JSON_CONTENT,
                            status=False,
                            detail=f"Validation failed for key '{key}' in object: {obj}",
                        )
            return create_response(
                name=ResponseName.VALID_EXAMPLE_JSON_CONTENT,
                status=True,
                detail=f"Validation failed for key '{key}' in object: {obj}",
            )
        else:
            example_file, nsample = self.get_flag_values(flags)
            if not Path(example_file).exists():
                return False, f"Example file not found: {example_file}"

            response = self._validate_csv_file(example_file, nsample)
            return response

    def _validate_csv_file(self, example_file, nsample):
        try:
            with open(
                example_file, mode="r", newline="", encoding="utf-8"
            ) as csvfile:
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
                if _row_len != nsample:
                    return create_response(
                        name=ResponseName.VALID_EXAMPLE_CSV_LENGTH,
                        status=False,
                        detail=f"Incorrect number of rows. Expected: {nsample}, Found: {len(row)}",
                    )
            return create_response(
                name=ResponseName.VALID_EXAMPLE_CSV_CONTENT,
                status=True,
            )
        except Exception as e:
            return False, f"Error reading CSV file: {str(e)}"

    def get_flag_values(self, args):
        try:
            file_flag_index = args.index("-f")
            nsamples_flag_index = args.index("-n")
            return args[file_flag_index + 1], args[nsamples_flag_index + 1]
        except (ValueError, IndexError):
            return None


def get_rule(rule_name, *args, **kwargs):
    rule_class = RULE_REGISTRY.get(rule_name)
    if not rule_class:
        raise ValueError(f"Rule '{rule_name}' is not registered.")
    return rule_class().check(*args, **kwargs)
