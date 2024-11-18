import csv
import json
import re
from datetime import datetime


def create_compound_input_csv(csv_path):
    input_data = [
        "COc1ccc2c(NC(=O)Nc3cccc(C(F)(F)F)n3)ccnc2c1",
        "O=C(O)c1ccccc1NC(=O)N1CCC(c2ccccc2C(F)(F)F)CC1",
        "O=C(O)c1cc(C(=O)O)c(C(=O)N(Cc2cccc(Oc3ccccc3)c2)[C@H]2CCCc3ccccc32)cc1C(=O)O",
        "Cc1ccc(N2CCN(Cc3nc4ccccc4[nH]3)CC2)cc1C",
        "Cc1cccc(NC(=O)CN2CCC(c3ccccn3)CC2)c1",
        "Clc1cccc(-c2nnnn2Cc2cccnc2)c1Cl",
        "CNC(=O)Nc1ccc2c(c1)CC[C@@]21OC(=O)N(CC(=O)N(Cc2ccc(F)cc2)[C@@H](C)C(F)(F)F)C1=O",
        "Cc1[nH]nc2ccc(-c3cncc(OC[C@@H](N)Cc4ccccc4)c3)cc12",
        "NCCCCCCCCCCNS(=O)(=O)c1cccc2c(Cl)cccc12",
    ]

    with open(csv_path, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Input"])
        for line in input_data:
            writer.writerow([line])


def save_as_json(result, output_file, remove_list=None):
    try:
        if remove_list:
            for item in remove_list:
                result = result.replace(item, "")

        stripped, formatted = re.split(r"\}\s*\{", result.strip()), []

        for i, line in enumerate(stripped):
            if i == 0:
                line = line + "}"
            elif i == len(stripped) - 1:
                line = "{" + line
            else:
                line = "{" + line + "}"
            formatted.append(line)

        _data = []
        for obj in formatted:
            try:
                _data.append(json.loads(obj))
            except json.JSONDecodeError as e:
                print(f"Skipping invalid JSON object: {obj}. Error: {e}")
                continue

        with open(output_file, "w") as f:
            json.dump(_data, f, indent=4)

    except Exception as e:
        raise ValueError(f"Error processing result: {e}")


def get_commands(model_id, config):
    return {
        "fetch": ["ersilia", "-v", "fetch", model_id]
        + (
            config.get("fetch_flags", "").split()
            if config.get("fetch_flags")
            else []
        ),
        "serve": ["ersilia", "-v", "serve", model_id]
        + (
            config.get("serve_flags", "").split()
            if config.get("serve_flags")
            else []
        ),
        "run": [
            "ersilia",
            "run",
            "-i",
            config["input_file"],
        ]
        + (
            ["-o", config["output_file"]]
            if config.get("output_file")
            and not config.get("output_redirection")
            else []
        )
        + (
            config.get("run_flags", "").split()
            if config.get("run_flags")
            else []
        ),
        "close": ["ersilia", "close"],
    }


def get_command_names(model_id, cli_type, config):
    return (
        list(get_commands(model_id, config).keys())
        if cli_type == "all"
        else [cli_type]
    )


def handle_error_logging(command, description, result, config, checkups=None):
    if config.get("log_error", False):
        log_path = f"{config.get('log_path')}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        with open(log_path, "w") as file:
            file.write(f"Command: {' '.join(command)}\n")
            file.write(f"Description: {description}\n")
            file.write(f"Error: {result}\n")
            if checkups:
                for check in checkups:
                    file.write(f"Check '{check['name']}': {check['status']}\n")
