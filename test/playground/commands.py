import psutil
import platform
import pytest
import subprocess
import time
import yaml
from pathlib import Path
from rich.text import Text

from .rules import get_rule
from .shared import results
from .utils import (
    get_command_names,
    get_commands,
    handle_error_logging,
)

config = yaml.safe_load(Path("config.yml").read_text())
delete_model = config.get("delete_model", False)
activate_docker = config.get("activate_docker", False)
runner = config.get("runner", "single")
cli_type = config.get("cli_type", "all")
output_file = config.get("output_files")
input_file = config.get("input_file")

if runner == "single":
    model_ids = [config["model_id"]]
else:
    model_ids = config["model_ids"]

base_path = Path.home() / "eos"
max_runtime_minutes = config.get("max_runtime_minutes", None)
from_github = "--from_github" in config.get("fetch_flags", "")
from_dockerhub = "--from_dockerhub" in config.get("fetch_flags", "")


def manage_docker(config):
    docker_activated = False

    if config and config.get("activate_docker"):
        system_platform = platform.system()

        if system_platform == "Linux":
            docker_status = subprocess.run(
                ["systemctl", "is-active", "--quiet", "docker"]
            )
            if docker_status.returncode != 0:
                subprocess.run(["systemctl", "start", "docker"], check=True)
        elif system_platform == "Darwin":  # macOS
            docker_status = subprocess.run(
                ["docker", "info"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
            if docker_status.returncode != 0:
                print(
                    "Docker is not running. Please start the Docker Desktop application."
                )
                raise RuntimeError(
                    "Docker is not running on macOS. Start Docker Desktop."
                )
        else:
            raise OSError(f"Unsupported platform: {system_platform}")

        docker_activated = True
    else:
        if platform.system() == "Linux":
            subprocess.run(["systemctl", "stop", "docker"], check=True)
        elif platform.system() == "Darwin":
            print("Stopping Docker programmatically is not supported on macOS.")
            raise RuntimeError("Cannot stop Docker programmatically on macOS.")

    return docker_activated


def execute_command(command, description="", dest_path=None, repo_path=None):
    docker_activated = manage_docker(config)
    (
        start_time,
        max_memory,
        success,
        result,
        checkups,
    ) = (
        time.time(),
        0,
        False,
        "",
        [],
    )

    proc = psutil.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    try:
        while proc.poll() is None:
            max_memory = max(max_memory, proc.memory_info().rss / (1024 * 1024))
            time.sleep(0.1)

        success = proc.returncode == 0
        stdout, stderr = proc.communicate()

        if success:
            result = stdout.decode()
        else:
            result = stderr.decode()

    except Exception as e:
        proc.kill()
        result = str(e)

        if config.get("log_error", False):
            handle_error_logging(command, description, result, config)

        pytest.fail(f"{description} '{' '.join(command)}' failed with error: {result}")

    checkups = apply_rules(command, description, dest_path, repo_path, config)

    rules_success = all(check["status"] for check in checkups)
    overall_success = success and rules_success

    if not overall_success and config.get("log_error", False):
        handle_error_logging(command, description, result, config, checkups)

    status_text = (
        Text("PASSED", style="green")
        if overall_success
        else Text("FAILED", style="red")
    )

    results.append(
        {
            "command": " ".join(command),
            "description": description,
            "time_taken": f"{(time.time() - start_time) / 60:.2f} min",
            "max_memory": f"{max_memory:.2f} MB",
            "status": status_text,
            "checkups": checkups,
            "activate_docker": docker_activated,
            "runner": config.get("runner"),
            "cli_type": config.get("cli_type"),
        }
    )

    return overall_success, (time.time() - start_time) / 60


def apply_rules(command, description, dest_path, repo_path, config):
    checkups = []
    try:
        if description == "fetch":
            if from_github:
                checkups.append(
                    get_rule(
                        "folder_exists",
                        folder_path=repo_path,
                        expected_status=True,
                    )
                )
            if from_dockerhub:
                checkups.append(
                    get_rule(
                        "folder_exists",
                        folder_path=dest_path,
                        expected_status=True,
                    )
                )
                checkups.append(
                    get_rule(
                        "dockerhub_status",
                        dest_path=dest_path,
                        expected_status=True,
                    )
                )
        elif "run" in description:
            file = next(
                (
                    t
                    for t in config["output_files"]
                    if Path(t).suffix.replace(".", "") in description
                ),
                None,
            )
            checkups.append(
                get_rule(
                    "file_exists",
                    file_path=file,
                    expected_status=True,
                )
            )

            inp_type = next(
                (t for t in config["input_types"] if t in description), None
            )

            checkups.append(
                get_rule(
                    "file_content_check",
                    file_path=file,
                    expected_input=inp_type,
                )
            )
        elif description == "serve":
            if from_dockerhub:
                checkups.append(
                    get_rule(
                        "dockerhub_status",
                        dest_path=dest_path,
                        expected_status=True,
                    )
                )
    except Exception as rule_error:
        handle_error_logging(command, description, rule_error, config, checkups)
        pytest.fail(f"Rule exception occurred: {rule_error}")

    return checkups


@pytest.fixture(scope="module", autouse=True)
def delete_model_command():
    if delete_model:
        for model_id in model_ids:
            dest_path = base_path / "dest" / model_id
            repo_path = base_path / "repository" / model_id
            delete_command = ["ersilia", "-v", "delete", model_id]
            success, _ = execute_command(
                delete_command,
                description="delete",
                dest_path=dest_path,
                repo_path=repo_path,
            )
            assert success, f"Delete command failed for model ID {model_id}"
            assert (
                not dest_path.exists()
            ), f"Destination folder for {model_id} still exists after delete"


@pytest.mark.parametrize(
    "model_id", model_ids if runner == "multiple" else [config.get("model_id")]
)
@pytest.mark.parametrize(
    "command_name", get_command_names(model_ids[0], cli_type, config)
)
def test_command(model_id, command_name):
    commands = get_commands(model_id, config)
    command = commands[command_name]
    dest_path = base_path / "dest" / model_id
    repo_path = base_path / "repository" / model_id

    success, time_taken = execute_command(
        command,
        description=command_name,
        dest_path=dest_path,
        repo_path=repo_path,
    )

    assert success, f"Command '{command_name}' failed for model ID {model_id}"

    if "run" in command_name and max_runtime_minutes is not None:
        assert (
            time_taken <= max_runtime_minutes
        ), f"Command '{command_name}' for model ID {model_id} exceeded max runtime of {max_runtime_minutes} minutes"
