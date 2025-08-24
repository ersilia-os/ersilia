import json
import os
import time
import traceback

import pytest
from click.testing import CliRunner
from rich.text import Text
from .shared import results
from .utils import (
    apply_rules,
    build_run_cmd,
    get_command,
    get_command_names,
    get_error,
    get_memory_usage,
    get_settings,
    handle_error_logging,
    handle_exception,
    manage_docker,
)

config = json.loads(os.getenv("CONFIG_DATA"))

(
    runner_mode, 
    cli, 
    silent, 
    models, 
    max_runtime_minutes, 
    base_path, 
    show_remark,
    model,
    host,
    activate_docker
) = get_settings(config)

def execute_command(command, description):
    def _execute(command, description):
        runner = CliRunner()
        docker_status = activate_docker
        if host == "local":
            docker_status = manage_docker(config)
        (
            start_time,
            max_memory,
            success,
            result,
            checkups,
            details
        ) = (
            time.time(),
            0,
            False,
            "",
            [],
            None
        )
        try:
            memory_before = get_memory_usage()
            res = runner.invoke(command[0], command[1])
            memory_after = get_memory_usage()

            max_memory = max(memory_before, memory_after)
            success = res.exit_code == 0
            result = res.output.strip()
            
            if not success:
                error_message = get_error(res)
                handle_exception(error_message, silent)
                handle_error_logging(command, description, error_message, config, checkups)
                pytest.fail(f"{description} '{command[2]}' failed:\n{error_message}")

        except Exception as e:
            handle_exception(e, silent)
            error_message = f"An exception occurred:\n{traceback.format_exc()}"
            handle_error_logging(command, description, error_message, config, checkups)
            pytest.fail(f"{description} '{command[2]}' failed with error:\n{error_message}")

        checkups = apply_rules(command, description, config, result)

        rules_success = all(check["status"] for check in checkups)
        overall_success = success and rules_success

        if not overall_success:
            details = [entry["detail"] for entry in checkups if "status" in entry]
            handle_error_logging(command, description, result, config, checkups)
        

        status_text = (
            Text("PASSED", style="green bold")
            if overall_success
            else Text("FAILED", style="red bold")
        )
        exec_time = (time.time() - start_time) / 60

        results.append(
            {
                "command": command[2],
                "time_taken": f"{exec_time:.2f} min",
                "max_memory": f"{max_memory:.2f} MB",
                "status": status_text,
                "checkups": checkups,
                "activate_docker": docker_status,
                "runner": runner_mode,
                "cli": cli,
                "show_remark": show_remark,
                "remark": result
            }
        )

        return overall_success, exec_time, details

    res = _execute(command, description)
    return res


@pytest.mark.parametrize(
    "model", models if runner_mode == "multiple" else [model]
)
@pytest.mark.parametrize("command_name", get_command_names(models[0], cli, config))
def test_command(model, command_name):
    commands = get_command(model, config)
    command = commands[command_name]
    if "run" in command_name:
        run_commands = build_run_cmd(config) 
        for _command in run_commands:
            success, time_taken, details = execute_command(
                _command,
                description=command_name,
            )
    else:
        success, time_taken, details = execute_command(
            command,
            description=command_name,
        )

    assert success, f"Command '{command_name}' failed for model ID {model} due to: {details}"

    if "run" in command_name and max_runtime_minutes is not None:
        assert time_taken <= max_runtime_minutes, (
            f"Command '{command_name}' for model ID {model} exceeded max \
                runtime of {max_runtime_minutes} minutes"
        )

    if "example" in command_name:
        if os.path.exists("tmp_example.csv"):
            os.remove("tmp_example.csv")
