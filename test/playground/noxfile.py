import json
import nox
import os
import shutil
import yaml
from pathlib import Path
from ersilia.default import EOS_PLAYGROUND

if not os.path.exists(EOS_PLAYGROUND):
    os.makedirs(EOS_PLAYGROUND)

NOX_CWD = Path.cwd()
NOX_PWD = Path(__file__).parent.parent.parent
CONFIG_PATH = NOX_CWD / "config.yml"
DEFAULT_CONFIG = yaml.safe_load(CONFIG_PATH.read_text())
PYTHON_VERSIONS = DEFAULT_CONFIG["settings"]["python_version"]
DEFAULT_BACKEND = DEFAULT_CONFIG["runtime"]["backend"]

nox.options.envdir = str(Path(EOS_PLAYGROUND) / ".nox")

test_packages = [
    "pytest", 
    "pytest-asyncio", 
    "pytest-benchmark", 
    "nox", 
    "rich",
    "fuzzywuzzy", 
    "scipy"
]

fetch_flags = (
    "from_github", 
    "from_dockerhub", 
    "from_s3", 
    "version", 
    "from_dir"
)

flagged_keys = { 
        "fetch": {*fetch_flags},
        "test": {*fetch_flags},
        "test": {*fetch_flags, "shallow", "deep", "as_json", "version"},
        "delete": {"all"},
        "catalog": {"more", "as-json", "local", "hub", "local"},
        "example": {"simple", "random", "predefined", "complete", "file_name", "n_samples"},
        "serve": {"enable-local-cache", "disable-local-cache"},
}

def parse_yaml(file_path):
    with open(file_path, "r") as f:
        return yaml.safe_load(f)


def preprocess_values(config_dict, flagged_keys):
    for key, known_values in flagged_keys.items():
        if key in config_dict:
            values = config_dict[key]
            if isinstance(values, str):  
                values = [values]
            if isinstance(values, list):
                config_dict[key] = [
                    f"--{val}" if val in known_values else val for val in values
                ]
    return config_dict

def parse_posargs(posargs, config_keys):
    parsed_dict = {}
    i = 0
    while i < len(posargs):
        arg = posargs[i]
        if arg.startswith("--"):
            key = arg[2:] 
            if key in config_keys:
                values = []
                i += 1
                while i < len(posargs) and not posargs[i].startswith("--"):
                    values.append(posargs[i])
                    i += 1
                parsed_dict[key] = values
            else:
                i += 1
        else:
            i += 1

    for key, value in parsed_dict.items():
        if isinstance(value, list):
            if len(value) == 1:
                if value[0].lower() in {"true", "false"}:
                    parsed_dict[key] = value[0].lower() == "true"
                elif value[0].isdigit():
                    parsed_dict[key] = int(value[0])
                else:
                    parsed_dict[key] = value[0]

    parsed_dict = preprocess_values(parsed_dict, flagged_keys)

    return parsed_dict

def get_all_subkeys(data):
    keys = []
    for key, value in data.items():
        if isinstance(value, dict):
            keys.extend(get_all_subkeys(value)) 
        else:
            keys.append(key)  
    return keys

def map_back_to_structure(parsed_data, default_config):
    reconstructed = {}
    top_level_keys = [key for key in default_config.keys() if isinstance(default_config[key], dict)]
    for section in top_level_keys:
        reconstructed[section] = {}
        
        for key in default_config[section].keys():
            if key in parsed_data:
                reconstructed[section][key] = parsed_data[key]
    
    return reconstructed

def replace_configs(default_config, override_config):
    replaced_config = default_config.copy()
    for key, value in replaced_config.items():
        if isinstance(value, dict): 
            if key in override_config and isinstance(override_config[key], dict):
                replaced_config[key] = replace_configs(value, override_config[key])
        elif key in override_config and override_config[key] is not None:
            replaced_config[key] = override_config[key]
    
    return replaced_config

def setup(session):
    session.log(f"Installing ersilia from source: {NOX_PWD}")
    session.install("-e", str(NOX_PWD))
    session.install(*test_packages)
    session.env["TEST_ENV"] = "true"

def run(session):
    default_config = parse_yaml(CONFIG_PATH)
    keys = get_all_subkeys(default_config)
    override_config = parse_posargs(session.posargs, keys)
    override_config = map_back_to_structure(override_config, default_config)
    final_config = replace_configs(default_config, override_config)
    session.env["CONFIG_DATA"] = json.dumps(final_config)
    session.cd(NOX_CWD)
    if final_config["settings"].get("silent"):
        session.run("pytest", "commands.py", "--disable-warnings")
    else:
        session.run("pytest", "-s", "commands.py", "--disable-warnings")

@nox.session(
    venv_backend=DEFAULT_BACKEND, python=PYTHON_VERSIONS, reuse_venv=True
)
def execute(session):
    setup(session)
    run(session)

@nox.session
def clean(session):
    session.log(f"Cleaning up {EOS_PLAYGROUND} and Config Data")
    if "CONFIG_DATA" in os.environ:
        del os.environ["CONFIG_DATA"]
    shutil.rmtree(EOS_PLAYGROUND, ignore_errors=True)
    session.log("Cleanup completed!")