import json
import os
import platform
import psutil
import pytest
import random
import requests
import subprocess
import click
import time
import traceback
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
from .rules import get_rule

file_path = Path(EOS_PLAYGROUND) / 'files'

CMD_DEPENDENCY_MAP = {
    "fetch": ["fetch"],
    "serve": ["fetch", "serve"],
    "run": ["fetch", "serve", "run"],
    "close": ["close"],
    "catalog": ["catalog"],
    "example": ["example"],
    "delete": ["delete"],
    "test": ["test"],
}

RAW_INPUT_URL = "https://github.com/ersilia-os/ersilia-model-hub-maintained-inputs/blob/main/compound/single/inp-000.csv"


def build_run_cmd(config):
    run_cli = []
    input_types = config["runtime"].get("input_types")
    input_types = input_types if isinstance(input_types, list) else [input_types]
    output_files = config["files"].get("outputs")
    output_files = output_files if isinstance(output_files, list) else [output_files]
    for _, inp_type in enumerate(input_types):
        for _, out_file in enumerate(output_files):
            _flag = ["-i", inp_type, "-o", out_file]
            out_file = os.path.join(file_path, out_file)
            inputs = get_inputs(config, inp_type)
            flags = ["-i", inputs, "-o", out_file]
            flag_description = f" | flag: {_flag}"
            run_cli.append((run_cmd(), flags, f"run_cmd(): Model id not required {flag_description}"))
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
        "github.com", 
        "raw.githubusercontent.com"
    ).replace("/blob/", "/")

    raw_inp_path = os.path.join(file_path, filename)
    num_samples  = int(config["settings"].get("num_samples"))
    input_file = os.path.join(
        file_path, 
        config["files"].get("input_single")
    )
    if not os.path.exists(raw_inp_path):
        response = requests.get(raw_url)
        if response.status_code == 200:
            os.makedirs(os.path.dirname(raw_inp_path), exist_ok=True)
            with open(raw_inp_path, "wb") as f:
                f.write(response.content)
        else:
            raise Exception(f"Failed to download file: {response.status_code}")

    with open(raw_inp_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    data = [line.strip().split(",") for line in lines[1:]]
    sampled_rows = random.sample(data, min(num_samples, len(data)))
    sampled_inputs = [row[1] for row in sampled_rows]

    with open(input_file, "w", encoding="utf-8") as f:
        f.write(lines[0])
        f.writelines(",".join(row) + "\n" for row in sampled_rows)

    return sampled_inputs, input_file


def get_commands_all(model_id, config):
    def build_command(cmd, extra_args=None, flag_key=None, model_id=None):
        extra_args = extra_args or []
        flags = config["flags"].get(flag_key, [])
        flag_description = f" | flag: {flags}" if flags else " | flag: Default"
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

    return {
        "fetch": build_command(fetch_cmd, flag_key="fetch", model_id=model_id),
        "serve": build_command(serve_cmd, model_id=model_id),
        "run": build_command(run_cmd, flag_key="run"),
        "close": build_command(close_cmd),
        "catalog": build_command(catalog_cmd, flag_key="catalog"),
        "example": build_command(example_cmd, flag_key="example", model_id=model_id),
        "delete": build_command(delete_cmd, flag_key="delete", model_id=model_id),
        "test": build_command(test_cmd, flag_key="test", model_id=model_id),
    }


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

    merged_dep = {dep for cmd in cli for dep in CMD_DEPENDENCY_MAP.get(cmd, [])}

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
    return (
        subprocess.run(
            ["docker", "info"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        ).returncode
        == 0
    )


def start_docker():
    system_platform = platform.system()

    if is_docker_running():
        print("Docker is already running. No action needed.")

    print("Starting Docker...")

    if system_platform == "Linux":
        subprocess.run(["sudo", "systemctl", "start", "docker"], check=True)

    elif system_platform == "Darwin":  # macOS
        print("Opening Docker Desktop...")
        subprocess.run(["open", "-a", "Docker"], check=True)

        print("Waiting for Docker to be fully operational...")
        for _ in range(30):
            if is_docker_running():
                print("Docker is ready!")
                return
            time.sleep(2)

        raise RuntimeError("Docker failed to start within the expected time.")

    else:
        raise OSError(f"Unsupported platform: {system_platform}")


def stop_docker():
    system_platform = platform.system()

    if not is_docker_running():
        print("Docker is already stopped. No action needed.")
        return

    print("Stopping Docker...")

    if system_platform == "Linux":
        subprocess.run(["sudo", "systemctl", "stop", "docker"], check=True)

    elif system_platform == "Darwin":  # macOS
        print("Stopping Docker programmatically is not supported on macOS.")
        raise RuntimeError("Cannot stop Docker programmatically on macOS.")

    else:
        raise OSError(f"Unsupported platform: {system_platform}")


def manage_docker(config):
    if not config or not config["settings"].get("activate_docker"):
        stop_docker()
        return False

    start_docker()
    return True


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


from datetime import datetime


def handle_error_logging(command, description, result, config, checkups):
    if config["settings"].get("log_error"):
        path = Path(EOS_PLAYGROUND) / 'logs'
        if not os.path.exists(path):
            os.makedirs(path)
        file_name = f"{description}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        log_path = os.path.join(path, file_name)
        with open(log_path, "w") as file:
            file.write(f"Command: {command}\n")
            file.write(f"Description: {description}\n")
            file.write(f"Error: {result}\n")
            if checkups:
                for check in checkups:
                    file.write(f"Check '{check['name']}': {check['status']} and Details: '{check['detail']}'\n")


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
    }

    try:
        if description in rule_mapping:
            checkups.extend(
                get_rule(
                    rule_mapping[description], 
                    command=command, 
                    config=config, 
                    std_out=std_out
            ))

    except Exception as rule_error:
        error_message = (
            f"Rule exception occurred: {rule_error}\n{traceback.format_exc()}"
        )
        handle_error_logging(
            command, 
            description, 
            rule_error, 
            config, 
            checkups
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
    model = model_config.get("single")

    return (
        runner_mode,
        cli,
        silent,
        models,
        max_runtime_minutes,
        base_path,
        show_remark,
        model,
    )

if __name__ == "__main__":
    import yaml
    config = yaml.safe_load(open("config.yml", "r"))
    commands = get_command("eos3b5e", config)
    fetch = commands["catalog"]
    print(fetch)