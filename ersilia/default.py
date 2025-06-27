import os
import shutil
from enum import Enum
from pathlib import Path

# ruff: noqa: D101, D102

# EOS environmental variables
EOS = os.path.join(str(Path.home()), "eos")
EOS_TMP = os.path.join(EOS, "temp")
if not os.path.exists(EOS):
    os.makedirs(EOS)
EOS_PLAYGROUND = os.path.join(EOS, "playground")
ROOT = os.path.dirname(os.path.realpath(__file__))
BENTOML_PATH = os.path.join(str(Path.home()), "bentoml")
CHECKSUM_NCHAR = 8
CLOUD_CACHE_CHUNK = 100_000
CONDA_ENV_YML_FILE = "environment.yml"
RUN_FILE = "run.sh"
DOCKERFILE_FILE = "Dockerfile"
GITHUB_ORG = "ersilia-os"
GITHUB_CONTENT_URL = f"https://raw.githubusercontent.com/{GITHUB_ORG}"
GITHUB_ERSILIA_REPO = "ersilia"
ERSILIA_MODEL_HUB_S3_BUCKET = "ersilia-model-hub"
ERSILIA_MODELS_S3_BUCKET = "ersilia-models"
ERSILIA_MODELS_ZIP_S3_BUCKET = "ersilia-models-zipped"
MODELS_JSON = "models_reannotated.json"
CONFIG_JSON = "config.json"
CREDENTIALS_JSON = "credentials.json"
INSTALL_STATUS_FILE = ".install.status"
DOCKER_BENTO_PATH = "/bento"
DOCKERHUB_ORG = "ersiliaos"
DOCKERHUB_LATEST_TAG = "latest"
DEFAULT_DOCKER_PLATFORM = "linux/amd64"
DEFAULT_MODEL_ID = "eos0zzz"
DEFAULT_VENV = "env"
DEFAULT_API_NAME = "run"
PACKMODE_FILE = "pack_mode.txt"
PACKMETHOD_FILE = "pack_method.txt"
CARD_FILE = "card.json"
UNPROCESSABLE_INPUT = "UNPROCESSABLE_INPUT"
DOTENV_FILE = ".env"
REDIS_DATA_VOLUME = "redis_data_volume"
API_SCHEMA_FILE = "api_schema.json"
MODEL_SIZE_FILE = "size.json"
DEFAULT_BATCH_SIZE = 100
DEFAULT_REDIS_MEMORY_USAGE_FRACTION = 0.3
FETCHED_MODELS_FILENAME = "fetched_models.txt"
MODEL_CONFIG_FILENAME = "config.json"
EXAMPLE_STANDARD_INPUT_CSV_FILENAME = "example_standard_input.csv"
EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME = "example_standard_output.csv"

PREDEFINED_EXAMPLE_INPUT_FILES = [
    "model/framework/examples/run_input.csv",
    "model/framework/examples/input.csv",
]

PREDEFINED_EXAMPLE_OUTPUT_FILES = [
    "model/framework/examples/output.csv",
    "model/framework/examples/run_output.csv",
]

PREDEFINED_COLUMN_FILE = "model/framework/columns/run_columns.csv"

PREDEFINED_EXAMPLE_FILES = (
    PREDEFINED_EXAMPLE_INPUT_FILES + PREDEFINED_EXAMPLE_OUTPUT_FILES
)
DEFAULT_ERSILIA_ERROR_EXIT_CODE = 1
METADATA_JSON_FILE = "metadata.json"
METADATA_YAML_FILE = "metadata.yml"
INSTALL_YAML_FILE = "install.yml"
SERVICE_CLASS_FILE = "service_class.txt"
MODEL_SOURCE_FILE = "model_source.txt"
APIS_LIST_FILE = "apis_list.txt"
INFORMATION_FILE = "information.json"
DOCKER_INFO_FILE = "from_dockerhub.json"
STATUS_JOSN = "status.json"
IS_FETCHED_FROM_HOSTED_FILE = "from_hosted.json"
DEFAULT_UDOCKER_USERNAME = "udockerusername"
DEFAULT_UDOCKER_PASSWORD = "udockerpassword"
DEFAULT_DOCKER_NETWORK_NAME = "ersilia_network"
DEFAULT_DOCKER_NETWORK_BRIDGE = "bridge"
# ERSILIA_RUNS_FOLDER = "ersilia_runs"
ALLOWED_API_NAMES = ["run", "train"]  # This can grow in the future based on needs
PACK_METHOD_FASTAPI = "fastapi"
PACK_METHOD_BENTOML = "bentoml"
BENTOML_APPROVED_PYTHON_VERSIONS = [
    "py36",
    "py37",
    "py38",
    "py39",
    "py310",
    "py311",
]
FASTAPI_APPROVED_PYTHON_VERSIONS = [
    "py38",
    "py39",
    "py310",
    "py311",
    "py312",
]
# Session and logging
SESSIONS_DIR = os.path.join(EOS, "sessions")
if not os.path.exists(SESSIONS_DIR):
    os.makedirs(SESSIONS_DIR, exist_ok=True)
SESSION_HISTORY_FILE = "history.txt"
SESSION_JSON = "session.json"
LOGS_DIR = "logs"
CONTAINER_LOGS_TMP_DIR = "_logs/tmp"
CONTAINER_LOGS_EOS_DIR = "_logs/eos"  # This is not used
LOGGING_FILE = "console.log"
CURRENT_LOGGING_FILE = "current.log"
SILENCE_FILE = ".silence.json"
VERBOSE_FILE = ".verbose.json"

H5_EXTENSION = ".h5"
H5_DATA_FILE = "data.h5"
REDIS_SERVER = "redis://redis:6379"
# Other
FEATURE_MERGE_PATTERN = "---"
HEADER_INDICATORS = [
    "key",
    "input",
    "smiles",
    "compound",
]  # TODO This list can be expanded. It is simply to detect headers easily.

# URLS
REDIS_EXPIRATION = 3600 * 24 * 7
REDIS_PORT = 6379
REDIS_CONTAINER_NAME = "redis"
REDIS_IMAGE = "redis:latest"
REDIS_HOST = "localhost"
ERSILIA_WEB_URL = "https://ersilia.io"
ERSILIA_MODEL_HUB_URL = "https://ersilia.io/model-hub"
AIRTABLE_MODEL_HUB_VIEW_URL = "https://airtable.com/shrNc3sTtTA3QeEZu"
S3_BUCKET_URL = "https://ersilia-models.s3.eu-central-1.amazonaws.com"
S3_BUCKET_URL_ZIP = "https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com"
INFERENCE_STORE_API_URL = (
    "https://5x2fkcjtei.execute-api.eu-central-1.amazonaws.com/dev/precalculations"
)
API_BASE = "https://hov95ejni7.execute-api.eu-central-1.amazonaws.com/dev/predict"
# EOS conda
_resolve_script = "conda_env_resolve.py"
resolve_script = os.path.join(EOS, _resolve_script)
if not os.path.exists(resolve_script):
    shutil.copyfile(
        os.path.join(ROOT, "utils", "supp", _resolve_script), resolve_script
    )


# Catalog table border constants
class TableConstants(str, Enum):
    TOP_LEFT = "┌"
    TOP_MIDDLE = "┬"
    TOP_RIGHT = "┐"
    HORIZONTAL = "─"
    VERTICAL = "│"
    MIDDLE_LEFT = "├"
    MIDDLE_MIDDLE = "┼"
    MIDDLE_RIGHT = "┤"
    BOTTOM_LEFT = "└"
    BOTTOM_MIDDLE = "┴"
    BOTTOM_RIGHT = "┘"
    CELL_PADDING = " "
    COLUMN_SEPARATOR = " | "


snippet = (
    """
# >>> ersilia >>>
# !! Contents within this block are managed by 'ersilia' !!
eosconda() {
    EOS_MODEL_ENV=$(python %s $2);
    conda $1 $EOS_MODEL_ENV
}

ersilia() {
    if [[ $1 == "conda" ]]; then
        eosconda "${@: 2}"
    elif [[ $1 == "auth" ]]; then
        gh auth "${@: 2}"
    else
        command ersilia "$@"
    fi
}
# <<< ersilia <<<
"""
    % resolve_script
)


def bashrc_path():
    home_path = Path.home()
    rc = os.path.join(home_path, ".bashrc")
    if os.path.exists(rc):
        return rc
    pr = os.path.join(home_path, ".bash_profile")
    if os.path.exists(pr):
        return pr


def has_profile_snippet():
    fn = bashrc_path()
    if not os.path.exists(fn):
        return False
    with open(fn, "r") as f:
        text = f.read()
    if snippet in text:
        return True
    else:
        return False


def bashrc_cli_snippet(overwrite=True):
    """Write a conda snippet in the user profile.

    This function writes on the user profile to create an executable to work
    with conda environments based on model identifiers.

    Motivation behind this function is to define an ersilia CLI.

    Args:
        - overwrite (bool): Overwrite the current bash profile file if the eosconda string is found.
    """
    fn = bashrc_path()
    if fn is None:
        return
    with open(fn, "r") as f:
        text = f.read()
    if snippet in text:
        if overwrite:
            text = text.split(snippet)[0] + text.split(snippet)[1]
        else:
            return
    with open(fn, "w") as f:
        f.write(text)
    with open(fn, "a+") as f:
        f.write(snippet)
