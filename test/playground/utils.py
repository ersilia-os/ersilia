import json
import os
import yaml
from datetime import datetime
from ersilia.io.input import ExampleGenerator
from pathlib import Path

config_path = Path("config.yml")
config = yaml.safe_load(config_path.read_text())

def get_commands(model_id, config):
    run = construct_run_cmd(config)
    return {
        "fetch": ["ersilia", "-v", "fetch", model_id]
        + (config.get("fetch_flags", "").split() if config.get("fetch_flags") else []),
        "serve": ["ersilia", "-v", "serve", model_id]
        + (config.get("serve_flags", "").split() if config.get("serve_flags") else []),
        **run,
        "close": ["ersilia", "close"],
    }


def get_command_names(model_id, cli_type, config):
    return (
        list(get_commands(model_id, config).keys()) if cli_type == "all" else [cli_type]
    )

def construct_run_cmd(config):
    run = {}
    for inp_type in config["input_types"]:
        for i, output_file in enumerate(config["output_files"]):
            out_type = os.path.splitext(output_file)[1].replace(".", "")
            inp_data = get_inputs(config, inp_type)
            key = f"run-{inp_type}-{out_type}"
            val = ["ersilia", "-v", "run", "-i", inp_data, "-o", output_file]
            run[key] = val
    return run


def get_inputs(config, types):
    example = ExampleGenerator(model_id=config["model_id"])
    samples = example.example(
        n_samples=config["number_of_input_samples"],
        file_name=None,
        simple=True,
        try_predefined=True,
    )

    samples = [sample["input"] for sample in samples]
    if types == "str":
        return samples[0]
    if types == "list":
        return json.dumps(samples)
    if types == "csv":
        example.example(
            n_samples=config["number_of_input_samples"],
            file_name=config["input_file"],
            simple=True,
            try_predefined=True,
        )
        return config["input_file"]


def handle_error_logging(command, description, result, config, checkups=None):
    if config.get("log_error", False):
        log_path = (
            f"{config.get('log_path')}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        )
        with open(log_path, "w") as file:
            file.write(f"Command: {' '.join(command)}\n")
            file.write(f"Description: {description}\n")
            file.write(f"Error: {result}\n")
            if checkups:
                for check in checkups:
                    file.write(f"Check '{check['name']}': {check['status']}\n")



if __name__ == "__main__":
    import subprocess
    import sys

    docker_activated = False

    if config and config.get("activate_docker"):
        if sys.platform == "linux":  # macOS
            try:
                docker_status = subprocess.run(
                    ["docker", "info"],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    check=True,
                )
                docker_activated = True
                print("Docker is running.")
            except subprocess.CalledProcessError:
                raise RuntimeError(
                    "Docker is not running. Please start Docker manually on macOS."
                )
        else:  # Assume Linux
            docker_status = subprocess.run(
                ["systemctl", "is-active", "--quiet", "docker"]
            )
            if docker_status.returncode != 0:
                subprocess.run(["systemctl", "start", "docker"], check=True)
            docker_activated = True
    # else:
    #     if sys.platform != "darwin":  # Only stop Docker on Linux
    #         subprocess.run(["systemctl", "stop", "docker"], check=True)

    print(docker_activated)
