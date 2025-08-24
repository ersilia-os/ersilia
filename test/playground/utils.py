import json
import os
import platform
import psutil
import shutil
import pytest
import requests
import subprocess
import click
import traceback
import tempfile
from datetime import datetime
from pathlib import Path
from ersilia.default import EOS_PLAYGROUND
from ersilia.cli.commands.example import example_cmd
from ersilia.cli.commands.catalog import catalog_cmd
from ersilia.cli.commands.fetch import fetch_cmd
from ersilia.cli.commands.serve import serve_cmd
from ersilia.cli.commands.close import close_cmd
from ersilia.cli.commands.run import run_cmd
from ersilia.cli.commands.test import test_cmd
from ersilia.cli.commands.delete import delete_cmd
from ersilia.cli import echo
from .rules import get_rule


file_path = Path(EOS_PLAYGROUND) / "files"

CMD_DEPENDENCY_MAP = {
    "fetch": ["fetch"],
    "serve": ["fetch", "serve"],
    "run": ["fetch", "serve", "run"],
    "close": ["serve", "close"],
    "catalog": ["catalog"],
    "example": ["serve", "example"],
    "delete": ["delete"],
    "test": ["test"],
}

RAW_INPUT_URL = "https://github.com/ersilia-os/ersilia-model-hub-maintained-inputs/blob/main/compound/single/inp-000.csv"


def build_run_cmd(config):
    run_cli = []
    input_types = config["runtime"].get("input_types")
    input_types = (
        input_types if isinstance(input_types, list) else [input_types]
    )
    output_files = config["files"].get("outputs")
    output_files = (
        output_files if isinstance(output_files, list) else [output_files]
    )
    for _, inp_type in enumerate(input_types):
        for _, out_file in enumerate(output_files):
            _flag = ["-i", inp_type, "-o", out_file]
            out_file = os.path.join(file_path, out_file)
            inputs = get_inputs(config, inp_type)
            flags = ["-i", inputs, "-o", out_file]
            flag_description = f" | flag: {_flag}"
            run_cli.append(
                (
                    run_cmd(),
                    flags,
                    f"run_cmd(): {flag_description}",
                )
            )
    return run_cli


def get_inputs(config, input_type):
    input_list, input_path = get_random_samples(config=config)
    if input_type == "str":
        return input_list[0]
    if input_type == "list":
        return json.dumps(input_list)
    if input_type == "csv":
        return input_path


def get_random_samples(config, filename="inp-000.csv"):
    raw_url = RAW_INPUT_URL.replace(
        "github.com", "raw.githubusercontent.com"
    ).replace("/blob/", "/")

    raw_inp_path = os.path.join(file_path, filename)
    num_samples = int(config["settings"].get("num_samples"))
    input_file = os.path.join(
        file_path, config["files"].get("input_single")
    )
    if not os.path.exists(raw_inp_path):
        response = requests.get(raw_url)
        if response.status_code == 200:
            os.makedirs(os.path.dirname(raw_inp_path), exist_ok=True)
            with open(raw_inp_path, "wb") as f:
                f.write(response.content)
        else:
            raise Exception(
                f"Failed to download file: {response.status_code}"
            )

    with open(raw_inp_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    data = [line.strip().split(",") for line in lines[1:]]
    input_list = [row[1] for row in data if len(row) > 1]
    input_list = input_list[500:]
    
    valid_input = [inp for inp in input_list]
    sampled = valid_input[:num_samples]

    with open(input_file, "w", encoding="utf-8") as f:
        f.write("input\n")
        f.writelines(row + "\n" for row in sampled)

    return sampled, input_file


def _set_n_sample_if_none(flags, flag_key):
    if (
        flags
        and "example" == flag_key
        and not any(flag in flags for flag in ("-n", "--n_samples"))
    ):
        flags = ["-n", 10] + flags
    return flags


def get_commands_all(model_id, config):
    def build_command(cmd, extra_args=None, flag_key=None, model_id=None):
        extra_args = extra_args or []
        flags = config["flags"].get(flag_key, [])
        flags = _set_n_sample_if_none(flags, flag_key)
        flag_description = (
            f" | flag: {flags}" if flags else " | flag: Default"
        )
        if model_id:
            return (
                cmd(),
                [model_id] + extra_args + flags,
                f"{cmd.__name__}: {model_id}{flag_description}",
            )
        return (
            cmd(),
            extra_args + flags,
            f"{cmd.__name__}: Model id not required{flag_description}",
        )

    tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
    tmp_file = os.path.join(tmp_folder, "example.csv")

    data = {
        "fetch": build_command(
            fetch_cmd, flag_key="fetch", model_id=model_id
        ),
        "serve": build_command(serve_cmd, flag_key="serve", model_id=model_id),
        "run": build_command(run_cmd, flag_key="run"),
        "catalog": build_command(catalog_cmd, flag_key="catalog"),
        "example": build_command(example_cmd, extra_args={"file_name": tmp_file}, flag_key="example"),
        "close": build_command(close_cmd),
        "delete": build_command(
            delete_cmd, flag_key="delete", model_id=model_id
        ),
        "test": build_command(test_cmd, flag_key="test", model_id=model_id),
    }

    shutil.rmtree(tmp_folder)

    return data


def extract_commands(keys, model_id, config):
    all_commands = get_commands_all(model_id, config)
    return {key: all_commands[key] for key in keys if key in all_commands}


def get_command(model_id, config):
    cli = config["settings"].get("cli")
    if not cli:
        raise ValueError("cli cannot be empty in the configuration.")

    if "all" in cli:
        return get_commands_all(model_id, config)

    if isinstance(cli, str):
        command = cli
        if command not in CMD_DEPENDENCY_MAP:
            raise ValueError(f"Unknown command: {command}")
        cmds = CMD_DEPENDENCY_MAP[command]
        return extract_commands(cmds, model_id, config)

    merged_dep = {
        dep for cmd in cli for dep in CMD_DEPENDENCY_MAP.get(cmd, [])
    }

    return extract_commands(merged_dep, model_id, config)


def get_command_names(model_id, cli, config):
    if "all" in cli:
        return list(get_commands_all(model_id, config).keys())
    if isinstance(cli, str):
        return CMD_DEPENDENCY_MAP[cli]
    if not cli:
        raise ValueError("cli cannot be empty in the configuration.")
    return cli


def is_docker_running():
    try:
        result = subprocess.run(
            ["docker", "info"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return result.returncode == 0
    except FileNotFoundError:
        echo("Docker command not found. Ensure Docker is installed and in the PATH.")
        return False

def start_docker():
    if is_docker_running():
        echo("Docker is already running. No action needed.", fg="green", bold=True)
        return

    echo("Starting Docker...", fg="green", bold=True)

    system_platform = platform.system()

    if system_platform == "Linux":
        subprocess.run(["sudo", "systemctl", "start", "docker"], check=True)

    elif system_platform == "Darwin":  # macOS
        if not is_docker_running():
            echo("Docker is not running. Please start Docker Desktop manually.", fg="red", bold=True)
            raise RuntimeError("Docker cannot be started programmatically on macOS.")
        else:
            echo("Docker is already running.", fg="green", bold=True)

    else:
        raise OSError(f"Unsupported platform: {system_platform}")

def stop_docker():
    system_platform = platform.system()

    if not is_docker_running():
        echo("Docker is already stopped. No action needed.", fg="green", bold=True)
        return

    echo("Stopping Docker...", fg="green", bold=True)

    if system_platform == "Linux":
        subprocess.run(["sudo", "systemctl", "stop", "docker"], check=True)

    elif system_platform == "Darwin":  # macOS
        echo("Stopping Docker programmatically is not supported on macOS.")
        raise RuntimeError("Cannot stop Docker programmatically on macOS.")

    else:
        raise OSError(f"Unsupported platform: {system_platform}")


def manage_docker(config):
    if not config or not config["settings"].get("activate_docker"):
        echo("Docker activation not required.", fg="yellow", bold=True)
        return False

    try:
        start_docker()
        return True
    except FileNotFoundError:
        echo("Docker is not installed or not in the PATH. Skipping Docker-related tasks.", fg="yellow", bold=True)
        return False


def get_memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)


def get_error(res):
    error_message = f"Command failed with exit code {res.exit_code}\n"
    error_message += f"Output:\n{res.output}\n"
    if res.exc_info:
        error_message += "Traceback:\n"
        error_message += "".join(traceback.format_exception(*res.exc_info))
    else:
        error_message += "No traceback available.\n"
    return error_message


def handle_exception(exc, verbose):
    if verbose:
        click.secho("Detailed traceback:", fg="yellow")
        traceback.print_exc()
    else:
        click.secho(f"Error: {str(exc)}", fg="red")
        click.secho("Run with --verbose for more details.", fg="cyan")


def handle_error_logging(command, description, result, config, checkups):
    if config["settings"].get("log_error"):
        path = Path(EOS_PLAYGROUND) / "logs"
        if not os.path.exists(path):
            os.makedirs(path)
        file_name = (
            f"{description}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        )
        log_path = os.path.join(path, file_name)
        with open(log_path, "w") as file:
            file.write(f"Command: {command}\n")
            file.write(f"Description: {description}\n")
            file.write(f"Error: {result}\n")
            if checkups:
                for check in checkups:
                    file.write(
                        f"Check '{check['name']}': {check['status']} and Details: '{check['detail']}'\n"
                    )


def apply_rules(command, description, config, std_out):
    checkups = []
    rule_mapping = {
        "catalog": "catalog_rule",
        "example": "example_rule",
        "delete": "delete_rule",
        "fetch": "fetch_rule",
        "serve": "serve_rule",
        "close": "close_rule",
        "run": "run_rule",
        "test": "test_rule",
    }

    try:
        if description in rule_mapping:
            checkups.extend(
                get_rule(
                    rule_mapping[description],
                    command=command,
                    config=config,
                    std_out=std_out,
                )
            )

    except Exception as rule_error:
        error_message = f"Rule exception occurred: {rule_error}\n{traceback.format_exc()}"
        handle_error_logging(
            command, description, rule_error, config, checkups
        )
        pytest.fail(error_message)

    return checkups


def get_settings(config):
    setting_config = config["settings"]
    model_config = config["models"]
    runner_mode = config["runtime"].get("runner")
    cli = setting_config.get("cli")
    silent = bool(setting_config.get("silent"))

    if runner_mode == "single":
        models = [model_config["single"]]
    else:
        models = model_config["multiple"]

    max_runtime_minutes = setting_config.get("max_runtime_minutes")
    base_path = Path.home() / "eos"
    show_remark = setting_config.get("show_remark")
    host = config["runtime"].get("host")
    model = model_config.get("single")
    activate_docker = setting_config.get("activate_docker")

    return (
        runner_mode,
        cli,
        silent,
        models,
        max_runtime_minutes,
        base_path,
        show_remark,
        model,
        host,
        activate_docker,
    )
