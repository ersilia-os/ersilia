import json
import os
import re
import requests
import yaml
from pathlib import Path
from ersilia.utils.conda import SimpleConda
from ersilia.utils.docker import SimpleDocker
from ersilia.utils.logging import logger
from ersilia.default import EOS, DOCKERHUB_ORG, SESSIONS_DIR, SESSION_JSON
from ersilia.utils.session import get_session_dir
from ersilia.publish.test import CheckService, IOService, STATUS_CONFIGS

RULE_REGISTRY = {}

config_path = Path("config.yml")
config = yaml.safe_load(config_path.read_text())


def create_response(name, status, detail=""):
    return {
        "name": name,
        "status": status,
        "detail": detail,
    }


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
                name="Dest Folder Exist",
                status=False,
                detail=f"Folder does not exist at: {path}",
            )
        return create_response(name="Dest Folder Exist", status=True)

    def _check_dest_folder_content(self, command):
        _, model, _, _ = get_configs(command)
        path = os.path.join(self.dest, model)
        if os.path.exists(path) and os.listdir(path):
            return create_response(name="Dest folder has content", status=True)
        return create_response(
            name="Dest folder has content",
            status=False,
            detail=f"No content exists at: {path}",
        )

    def _check_docker_image(self, command):
        _, model, source, version = get_configs(command)
        version = version if version else "latest"
        if "dockerhub" in source:
            exist = self.sd.exists(org=DOCKERHUB_ORG, img=model, tag=version)
            return create_response(name="Docker Image Exist", status=exist)

    def _check_conda_env(self, command):
        _, model, source, _ = get_configs(command)
        if "github" in source or "s3" in source:
            exist = self.sc.exists(model)
            return create_response(name="Conda venv exists", status=exist)


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
            return create_response(name="Session required files exist", status=True)
        return create_response(
            name="Session required files exist",
            status=False,
            detail=f"Session required files do not exist at: {self.current_session_dir}",
        )

    def check_service_class(self):
        service_class = ("pulled_docker", "conda")
        session_file_path = os.path.join(self.current_session_dir, SESSION_JSON)

        if not os.path.exists(session_file_path):
            return create_response(
                name="Session file exists",
                status=False,
                detail=f"Session file does not exist at: {session_file_path}",
            )

        try:
            with open(session_file_path, "r") as file:
                data = json.load(file)
                correct_srv_class = data.get("service_class") in service_class
                return create_response(
                    name=f"Session service class is correct: {data.get('service_class')}",
                    status=correct_srv_class,
                )
        except (json.JSONDecodeError, KeyError):
            return create_response(
                name="Session service class is correct",
                status=False,
                detail=f"Session file is not valid at: {session_file_path}",
            )

    def extract_url(self, text):
        pattern = r"http://[a-zA-Z0-9.-]+:\d+"
        match = re.search(pattern, text)
        url = match.group(0)
        response = requests.get(url)
        status_ok = response.status_code == 200
        return create_response(name="Session url is valid", status=status_ok)


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
            return [create_response(name="file content is valid", status=True)]
        else:
            return [
                create_response(
                    name="file content is not valid", status=False, detail=res[1]
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
                name="Repo Folder Not Exist",
                status=False,
                detail=f"Folder does not exist at: {path}",
            )
        return create_response(name="Repo Folder Exist", status=True)

    def _check_dest_folder(self, command):
        _, model, _, _ = get_configs(command)
        path = os.path.join(self.dest, model)
        if os.path.exists(path):
            return create_response(
                name="Dest Folder Not Exist",
                status=False,
                detail=f"Folder does not exist at: {path}",
            )
        return create_response(name="Dest Folder Exist", status=True)

    def _check_containers_images_removed(self, command):
        _, model, source, version = get_configs(command)
        version = version if version else "latest"
        img_name = f"{DOCKERHUB_ORG}/{model}:{version}"
        containers = self.sd.containers(only_run=False)
        if containers:
            exists = img_name in containers.values()
            return create_response(name="Containers Removed", status=not exists)
        return create_response(name="Containers Removed", status=True)

    def _check_docker_image(self, command):
        _, model, source, version = get_configs(command)
        version = version if version else "latest"
        if "dockerhub" in source:
            exist = self.sd.exists(org=DOCKERHUB_ORG, img=model, tag=version)
            return create_response(name="Docker Image Not Exist", status=exist)
        return create_response(name="Docker Image Not Exist", status=True)

    def _check_conda_env(self, command):
        _, model, source, _ = get_configs(command)
        if self.sc.exists(model):
            return create_response(name="Conda venv not exists", status=False)
        return create_response(name="Conda venv not exists", status=True)


def get_rule(rule_name, *args, **kwargs):
    rule_class = RULE_REGISTRY.get(rule_name)
    if not rule_class:
        raise ValueError(f"Rule '{rule_name}' is not registered.")
    return rule_class().check(*args, **kwargs)
