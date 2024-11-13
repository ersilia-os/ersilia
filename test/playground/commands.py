import pytest
import subprocess
import time
import psutil
import yaml
from pathlib import Path
from .shared import results

config = yaml.safe_load(Path("config.yml").read_text())
delete_model = config.get("delete_model", False)
activate_docker = config.get("activate_docker", False)
runner = config.get("runner", "single")
cli_type = config.get("cli_type", "all")

if runner == "single":
    model_ids = [config["model_id"]]
else:
    model_ids = config["model_ids"]
base_path = Path.home() / "eos"


def get_commands(model_id):
    return {
        "fetch": ["ersilia", "-v", "fetch", model_id]
        + ([config["fetch_flags"]] if config.get("fetch_flags") else []),
        "serve": ["ersilia", "-v", "serve", model_id]
        + ([config["serve_flags"]] if config.get("serve_flags") else []),
        "run": ["ersilia", "run", "-i", config["input_file"]]
        + (["-o", config["output_file"]] if config.get("output_file") else [])
        + ([config["run_flags"]] if config.get("run_flags") else []),
        "close": ["ersilia", "close"],
    }


def get_command_names(model_id):
    commands = get_commands(model_id)
    return list(commands.keys()) if cli_type == "all" else [cli_type]


max_runtime_minutes = config.get("max_runtime_minutes", None)
from_github = "--from_github" in config.get("fetch_flags", "")
from_dockerhub = "--from_dockerhub" in config.get("fetch_flags", "")


def check_folder_exists_and_not_empty(path):
    return path.exists() and any(path.iterdir())


def check_dockerhub_status(expected_status, dest_path):
    dockerhub_file = dest_path / "from_dockerhub.json"
    if dockerhub_file.exists():
        with open(dockerhub_file, "r") as f:
            content = f.read()
        return f'"docker_hub": {str(expected_status).lower()}' in content
    return False


def execute_command(command, description="", dest_path=None, repo_path=None):
    docker_activated = False
    if activate_docker:
        docker_status = subprocess.run(
            ["systemctl", "is-active", "--quiet", "docker"]
        )
        if docker_status.returncode != 0:
            subprocess.run(["systemctl", "start", "docker"], check=True)
        docker_activated = True
    else:
        subprocess.run(["systemctl", "start", "docker"], check=True)


    start_time = time.time()
    proc = psutil.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    max_memory = 0
    success = False
    result = ""
    checkups = []

    try:
        while proc.poll() is None:
            max_memory = max(
                max_memory, proc.memory_info().rss / (1024 * 1024)
            )  # MB
            time.sleep(0.1)  # Poll every 100 ms
    except Exception as e:
        proc.kill()
        pytest.fail(
            f"{description} '{' '.join(command)}' failed with error: {e}"
        )

    time_taken = (time.time() - start_time) / 60
    success = proc.returncode == 0
    result = (
        proc.communicate()[0].decode()
        if success
        else proc.stderr.read().decode()
    )

    if description == "fetch":
        if from_github:
            checkups.append(
                {
                    "name": "Repo folder and files Exists",
                    "status": check_folder_exists_and_not_empty(repo_path),
                }
            )
            checkups.append(
                {
                    "name": "Pulled from docker",
                    "status": check_dockerhub_status(False, dest_path),
                }
            )
        if from_dockerhub:
            checkups.append(
                {
                    "name": "Dest folder exists",
                    "status": check_folder_exists_and_not_empty(dest_path),
                }
            )
            checkups.append(
                {
                    "name": "Pulled from docker",
                    "status": check_dockerhub_status(True, dest_path),
                }
            )
    elif description == "run":
        output_file = Path(config["output_file"])
        checkups.append(
            {
                "name": "Output File Not Nulls",
                "status": output_file.exists()
                and "null" not in output_file.read_text(),
            }
        )
    elif description == "serve":
        if from_github:
            checkups.append(
                {
                    "name": "Pulled from docker",
                    "status": check_dockerhub_status(False, dest_path),
                }
            )
        if from_dockerhub:
            checkups.append(
                {
                    "name": "Pulled from docker",
                    "status": check_dockerhub_status(True, dest_path),
                }
            )

    checkups.append({"name": "Docker Activated", "status": docker_activated})

    time_color = (
        "\033[91m"
        if (description == "run" and time_taken > max_runtime_minutes)
        else ""
    )
    error_message = (
        result if not success else "\033[92mExecuted normally\033[0m"
    )

    results.append(
        {
            "command": " ".join(command),
            "description": description,
            "time_taken": f"{time_taken:.2f} min",
            "max_memory": f"{max_memory:.2f} MB",
            "status": result if not success else "Executed normally",
            "checkups": checkups,
            "activate_docker": activate_docker,
            "runner": runner,
            "cli_type": cli_type,
        }
    )

    return success, time_taken


@pytest.fixture(scope="module", autouse=True)
def delete_model_command():
    """Deletes models if delete_model is True."""
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
@pytest.mark.parametrize("command_name", get_command_names(model_ids[0]))
def test_command(model_id, command_name):
    commands = get_commands(model_id)
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

    if command_name == "run" and max_runtime_minutes is not None:
        assert (
            time_taken <= max_runtime_minutes
        ), f"Command '{command_name}' for model ID {model_id} exceeded max runtime of {max_runtime_minutes} minutes"
