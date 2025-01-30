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
from ersilia.default import EOS, DOCKERHUB_ORG, SESSIONS_DIR, SESSION_JSON
from ersilia.utils.session import get_session_dir
from ersilia.publish.test import CheckService, IOService, STATUS_CONFIGS


# Enum for response names
class ResponseName(Enum):
    DEST_FOLDER_EXIST = "Dest Folder Exist"
    DEST_FOLDER_HAS_CONTENT = "Dest folder has content"
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
        raise NotImplementedError("Each rule must implement a check method.")


def register_rule(name):
    def decorator(cls):
        RULE_REGISTRY[name] = cls
        return cls
    return decorator


@register_rule("fetch_rule")
class FetchRule(CommandRule):
    def __init__(self):
        self.sc = SimpleConda()
        self.sd = SimpleDocker()
        self.dest = os.path.join(EOS, "dest")
        self.repos = os.path.join(EOS, "repository")
        self.checkups = []

    def check(self, command, config, std_out):
        self.checkups.append(self._check_folders(command))
        img_chks = self._check_docker_image(command)
        if img_chks:
            self.checkups.append(img_chks)
        self.checkups.append(self._check_dest_folder_content(command))
        env_chks = self._check_conda_env(command)
        if env_chks:
            self.checkups.append(env_chks)
        return self.checkups

    def _check_folders(self, command):
        _, model, _, _ = get_configs(command)
        path = os.path.join(self.dest, model)
        if not os.path.exists(path):
            return create_response(
                name=ResponseName.DEST_FOLDER_EXIST,
                status=False,
                detail=f"Folder does not exist at: {path}",
            )
        return create_response(name=ResponseName.DEST_FOLDER_EXIST, status=True)

    def _check_dest_folder_content(self, command):
        _, model, _, _ = get_configs(command)
        path = os.path.join(self.dest, model)
        if os.path.exists(path) and os.listdir(path):
            return create_response(
                name=ResponseName.DEST_FOLDER_HAS_CONTENT, status=True
            )
        return create_response(
            name=ResponseName.DEST_FOLDER_HAS_CONTENT,
            status=False,
            detail=f"No content exists at: {path}",
        )

    def _check_docker_image(self, command):
        _, model, source, version = get_configs(command)
        version = version if version else "latest"
        if "dockerhub" in source:
            exist = self.sd.exists(org=DOCKERHUB_ORG, img=model, tag=version)
            return create_response(name=ResponseName.DOCKER_IMAGE_EXIST, status=exist)

    def _check_conda_env(self, command):
        _, model, source, _ = get_configs(command)
        if "github" in source or "s3" in source:
            exist = self.sc.exists(model)
            return create_response(name=ResponseName.CONDA_ENV_EXISTS, status=exist)


@register_rule("serve_rule")
class ServeRule(CommandRule):
    def __init__(self):
        self.session_dir = SESSIONS_DIR
        self.checkups = []
        self.current_session_dir = get_session_dir()

    def check(self, command, config, std_out):
        print(f"Session dir: {get_session_dir()}")
        self.checkups.append(self.get_latest_folder_and_check_files(command))
        self.checkups.append(self.check_service_class())
        self.checkups.append(self.extract_url(std_out))
        return self.checkups

    def get_latest_folder_and_check_files(self, command):
        _, model, _, _ = get_configs(command)
        required_files = {
            "console.log",
            "current.log",
            "_logs",
            "logs",
            f"{model}.pid",
            SESSION_JSON,
        }

        existing_files = set(os.listdir(self.current_session_dir))
        required_files_exists = required_files.issubset(existing_files)
        if required_files_exists:
            return create_response(
                name=ResponseName.SESSION_REQUIRED_FILES_EXIST, status=True
            )
        return create_response(
            name=ResponseName.SESSION_REQUIRED_FILES_EXIST,
            status=False,
            detail=f"Session required files do not exist at: {self.current_session_dir}",
        )

    def check_service_class(self):
        service_class = ("pulled_docker", "conda")
        session_file_path = os.path.join(self.current_session_dir, SESSION_JSON)

        if not os.path.exists(session_file_path):
            return create_response(
                name=ResponseName.SESSION_FILE_EXISTS,
                status=False,
                detail=f"Session file does not exist at: {session_file_path}",
            )

        try:
            with open(session_file_path, "r") as file:
                data = json.load(file)
                correct_srv_class = data.get("service_class") in service_class
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
        inp_type = self._get_input_type(output_files)
        res = self.check_service.validate_file_content(output_files, inp_type)
        if res[-1] == str(STATUS_CONFIGS.PASSED):
            return [create_response(name=ResponseName.FILE_CONTENT_VALID, status=True)]
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
        _, model, _, _ = get_configs(command)
        path = os.path.join(self.repos, model)
        if os.path.exists(path):
            return create_response(
                name=ResponseName.REPO_FOLDER_NOT_EXIST,
                status=False,
                detail=f"Folder does not exist at: {path}",
            )
        return create_response(name=ResponseName.REPO_FOLDER_EXIST, status=True)

    def _check_dest_folder(self, command):
        _, model, _, _ = get_configs(command)
        path = os.path.join(self.dest, model)
        if os.path.exists(path):
            return create_response(
                name=ResponseName.DEST_FOLDER_NOT_EXIST,
                status=False,
                detail=f"Folder does not exist at: {path}",
            )
        return create_response(name=ResponseName.DEST_FOLDER_EXIST, status=True)

    def _check_containers_images_removed(self, command):
        _, model, source, version = get_configs(command)
        version = version if version else "latest"
        img_name = f"{DOCKERHUB_ORG}/{model}:{version}"
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
        pass

    def check(self, command, config, std_out):
        is_valid, detail = self.validate_json_data(std_out, command)
        return [
            create_response(
                name=ResponseName.CATALOG_JSON_CONTENT_VALID,
                status=is_valid,
                detail=detail,
            )
        ]

    def validate_json_data(self, std_out, command):
        flags, _, _, _ = get_configs(command)
        data = json.loads(std_out)
        if "--as-json" in flags:
            for obj in data:
                for key, value in obj.items():
                    if value is None or value == "" or value == [] or value == {}:
                        return (
                            False,
                            f"Validation failed for key '{key}' in object: {obj}",
                        )
            return True, ""


@register_rule("example_rule")
class ExampleRule(CommandRule):
    def __init__(self):
        pass

    def check(self, command, config, std_out):
        is_valid, detail = self.validate_json_data(std_out, command, config)
        return [
            create_response(
                name=ResponseName.VALID_EXAMPLE_JSON_CONTENT,
                status=is_valid,
                detail=detail,
            )
        ]

    def validate_json_data(self, std_out, command, config):
        flags, _, _, _ = get_configs(command)
        data = json.loads(std_out)
        if config["files"].get("example_file") not in flags:
            for obj in data:
                for key, value in obj.items():
                    if value is None or value == "" or value == [] or value == {}:
                        return (
                            False,
                            f"Validation failed for key '{key}' in object: {obj}",
                        )
            return True, ""


def get_rule(rule_name, *args, **kwargs):
    rule_class = RULE_REGISTRY.get(rule_name)
    if not rule_class:
        raise ValueError(f"Rule '{rule_name}' is not registered.")
    return rule_class().check(*args, **kwargs)
