import pytest
import subprocess
import time
import psutil
import yaml
from pathlib import Path
from .shared import results

config = yaml.safe_load(Path("config.yml").read_text())
delete_model = config.get("delete_model", False)

commands = {
    "fetch": ["ersilia", "-v", "fetch", config["model_id"]]
    + ([config["fetch_flags"]] if config.get("fetch_flags") else []),
    "serve": ["ersilia", "-v", "serve", config["model_id"]]
    + ([config["serve_flags"]] if config.get("serve_flags") else []),
    "run": ["ersilia", "run", "-i", config["input_file"]]
    + (["-o", config["output_file"]] if config.get("output_file") else [])
    + ([config["run_flags"]] if config.get("run_flags") else []),
    "close": ["ersilia", "close"],
}

max_runtime_minutes = config.get("max_runtime_minutes", None)

def execute_command(command, description=""):
    start_time = time.time()
    proc = psutil.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    max_memory = 0
    success = False
    result = ""

    try:
        while proc.poll() is None:
            max_memory = max(max_memory, proc.memory_info().rss / (1024 * 1024))  # MB
            time.sleep(0.1)  # Poll every 100 ms
    except Exception as e:
        proc.kill()
        pytest.fail(f"{description} '{' '.join(command)}' failed with error: {e}")

    time_taken = (time.time() - start_time) / 60
    success = proc.returncode == 0
    result = proc.communicate()[0].decode() if success else proc.stderr.read().decode()

    time_color = (
        "\033[91m"
        if (description == "run" and time_taken > max_runtime_minutes)
        else ""
    )
    error_message = result if not success else "\033[92mExecuted normally\033[0m"

    results.append(
        {
            "command": " ".join(command),
            "description": description,
            "time_taken": f"{time_color}{time_taken:.2f} min\033[0m",
            "max_memory": f"{max_memory:.2f} MB",
            "status": error_message,
        }
    )

    return success, time_taken


@pytest.fixture(scope="module", autouse=True)
def delete_model_command():
    if delete_model:
        delete_command = ["ersilia", "-v", "delete", config["model_id"]]
        success, _ = execute_command(delete_command, description="delete")
        assert success, "Delete command failed"


@pytest.mark.parametrize("command_name", commands.keys())
def test_command(command_name):
    command = commands[command_name]
    success, time_taken = execute_command(command, description=command_name)

    assert success, f"Command '{' '.join(command)}' failed"

    if command_name == "run" and max_runtime_minutes is not None:
        assert (
            time_taken <= max_runtime_minutes
        ), f"Command '{' '.join(command)}' exceeded max runtime of {max_runtime_minutes} minutes"