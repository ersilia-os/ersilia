import os
from dataclasses import dataclass
from enum import Enum
from typing import List

from ....default import (
    EOS_TMP,
    DOCKERFILE_FILE,
    INSTALL_YAML_FILE,
    METADATA_JSON_FILE,
    METADATA_YAML_FILE,
    EOS_TMP,
    PREDEFINED_EXAMPLE_INPUT_FILES,
    PREDEFINED_EXAMPLE_OUTPUT_FILES,
    PREDEFINED_COLUMN_FILE,
    RUN_FILE,
)


class Options(Enum):
    """
    Enum for different options.
    """

    NUM_SAMPLES = 5
    BASE = "base"
    DEEP_INPUT = "deep_input.csv"
    DEEP_OUTPUT = "deep_output.csv"
    OUTPUT_CSV = "result.csv"
    EXAMPLE_CSV = "example.csv"
    INPUT_CSV = "input.csv"
    OUTPUT1_CSV = "output1.csv"
    OUTPUT2_CSV = "output2.csv"
    OUTPUT_FILES = [
        "file.csv",
        "file.h5",
        "file.json",
    ]
    OUTPUT_FILES_TEST = [
        "file.csv",
    ]
    INPUT_TYPES = ["str", "list", "csv"]
    INPUT_TYPES_TEST = ["csv"]

    def __getattribute__(self, name):
        file_path = os.path.join(EOS_TMP, "files")
        if not os.path.exists(file_path):
            os.makedirs(file_path, exist_ok=True)
        value = super().__getattribute__(name)
        if (
            name == "value"
            and isinstance(value, str)
            and not value.startswith(EOS_TMP)
            and any(value.endswith(ext) for ext in (".csv", "h5", "json"))
        ):
            return os.path.join(file_path, value)
        elif (
            name == "value"
            and isinstance(value, list)
            and any(
                path.endswith(ext) for path in value for ext in (".csv", "h5", "json")
            )
        ):
            return [
                os.path.join(file_path, file) if not file.startswith(EOS_TMP) else file
                for file in value
            ]
        return value


class STATUS_CONFIGS(Enum):
    """
    Enum for status configurations.
    """

    PASSED = ("PASSED", "green", "✔")
    NOT_PRESENT = ("NOT PRESENT", "yellow", "✘")
    FAILED = ("FAILED", "red", "✘")
    WARNING = ("WARNING", "yellow", "⚠")
    SUCCESS = ("SUCCESS", "green", "★")
    NA = ("N/A", "dim", "~")

    def __init__(self, label, color, icon):
        self.label = label
        self.color = color
        self.icon = icon

    def __str__(self):
        return f"[{self.color}]{self.icon} {self.label}[/{self.color}]"


class Checks(Enum):
    """
    Enum for different check types.
    """

    MODEL_CONSISTENCY = "Check Consistency of Model Output"
    CSV_CONTENT_CHECK = "CSV File Content Check"
    IMAGE_SIZE = "Image Size Mb"
    PREDEFINED_EXAMPLE = "Check Predefined Example Input"
    ENV_SIZE = "Environment Size Mb"
    DIR_SIZE = "Directory Size Mb"
    # messages
    SIZE_CACL_SUCCESS = "Size Successfully Calculated"
    FETCH_FAILS = "Fetch Status"
    EMPTY_COLUMNS = "Empty Column Found"
    COLUMN_MISMATCH = "Column mismatch check"
    SIZE_CACL_FAILED = "Size Calculation Failed"
    INCONSISTENCY = "Inconsistent Output Detected"
    CONSISTENCY = "Model Output Was Consistent"
    RUN_BASH = "RMSE-MEAN"
    TOTAL_DIR_SIZE = "Total directory size"
    COLUMN_NAME_VALIDITY = "Columns"
    COLUMN_CHECK_SUCCESS = "Columns coincides with run_columns"
    COLUMN_CHECK_FAILURE = "Columns not coincide with run_columns"
    SIMPLE_MODEL_RUN = "Simple Model Run"
    SIMPLE_MODEL_RUN_COLUMNS = "Simple Model Run Columns"
    DEPENDENCY_PINNED = "Checking package versions and file structure"


class TableType(Enum):
    """
    Enum for different table types.
    """

    MODEL_INFORMATION_CHECKS = "Model Metadata Checks"
    FETCH_STATUS_SURFACE = "Model Fetching Check"
    MODEL_FILE_CHECKS = "Model File Checks"
    MODEL_DIRECTORY_SIZES = "Directory Size Check"
    MODEL_SIZES = "Model Size Check"
    RUNNER_CHECKUP_STATUS = "Runner Checkup Status"
    FINAL_RUN_SUMMARY = "Test Run Summary"
    DEPENDECY_COLUMN_CHECK = "Dependency and Column Value Checks"
    COMPUTATIONAL_PERFORMANCE = "Computational Performance Check"
    SHALLOW_CHECK_SUMMARY = "Model Output Consistency Check"
    CONSISTENCY_BASH = "Consistency Summary Between Ersilia and Bash Execution Outputs"
    MODEL_OUTPUT = "Input Output Check"

    INSPECT_SUMMARY = "Inspect Summary"
    MODEL_RUN_CHECK = "Model Run Check"


@dataclass
class TableConfig:
    """
    Configuration for a table.
    """

    title: str
    headers: List[str]


TABLE_CONFIGS = {
    TableType.MODEL_INFORMATION_CHECKS: TableConfig(
        title="\nMetadata Checks", headers=["Check", "Details", "Status"]
    ),
    TableType.MODEL_FILE_CHECKS: TableConfig(
        title="\nModel File Checks", headers=["Check", "Details", "Status"]
    ),
    TableType.MODEL_DIRECTORY_SIZES: TableConfig(
        title="\nDirectory Size Check", headers=["Check", "Details", "Size"]
    ),
    TableType.RUNNER_CHECKUP_STATUS: TableConfig(
        title="\nRunner Checkup Status",
        headers=["Runner", "Status"],
    ),
    TableType.FINAL_RUN_SUMMARY: TableConfig(
        title="\nTest Run Summary", headers=["Check", "Status"]
    ),
    TableType.DEPENDECY_COLUMN_CHECK: TableConfig(
        title="\nFile Validity Check", headers=["Check", "Details", "Status"]
    ),
    TableType.COMPUTATIONAL_PERFORMANCE: TableConfig(
        title="\nComputational Performance Summary", headers=["Check", "Status"]
    ),
    TableType.SHALLOW_CHECK_SUMMARY: TableConfig(
        title="\nModel Output Consistency Check",
        headers=["Check", "Details", "Status"],
    ),
    TableType.MODEL_OUTPUT: TableConfig(
        title="\nInput Output Check",
        headers=["Check", "Detail", "Status"],
    ),
    TableType.MODEL_SIZES: TableConfig(
        title="\nModel Size Check",
        headers=["Check", "Detail", "Status"],
    ),
    TableType.CONSISTENCY_BASH: TableConfig(
        title="\nConsistency Summary Between Ersilia and Bash Execution Outputs",
        headers=["Check", "Result", "Status"],
    ),
    TableType.INSPECT_SUMMARY: TableConfig(
        title="\nInspect Summary", headers=["Check", "Status"]
    ),
    TableType.MODEL_RUN_CHECK: TableConfig(
        title="\nModel Run Check", headers=["Check", "Details", "Status"]
    ),
    TableType.FETCH_STATUS_SURFACE: TableConfig(
        title="\nModel Fetching Check", headers=["Check", "Status"]
    )
}


RUN_FILE = f"model/framework/{RUN_FILE}"

BENTOML_FILES = [
    DOCKERFILE_FILE,
    METADATA_JSON_FILE,
    RUN_FILE,
    "src/service.py",
    "pack.py",
    "README.md",
    "LICENSE",
]

TIMEOUT_SECONDS = 30 * 60 

ERSILIAPACK_BACK_FILES = [
    DOCKERFILE_FILE,
    METADATA_JSON_FILE,
    RUN_FILE,
    PREDEFINED_EXAMPLE_INPUT_FILES[0],
    PREDEFINED_EXAMPLE_OUTPUT_FILES[-1],
    PREDEFINED_COLUMN_FILE,
    "src/service.py",
    "pack.py",
    "README.md",
    "LICENSE",
]

ERSILIAPACK_FILES = [
    INSTALL_YAML_FILE,
    METADATA_YAML_FILE,
    PREDEFINED_EXAMPLE_INPUT_FILES[0],
    PREDEFINED_EXAMPLE_OUTPUT_FILES[-1],
    PREDEFINED_COLUMN_FILE,
    RUN_FILE,
    "README.md",
    "LICENSE",
]


COMMON_FILES = [
    RUN_FILE,
    "README.md",
    "LICENSE",
]

BENTOML_FOLDERS = ["model", "src", ".github"]

ERSILIAPACK_FOLDERS = ["model", ".github"]

GIT_AND_DOCKER_IGNORE = [".dockerignore", ".gitignore", ".gitattributes"]

# Base URL for the Ersilia OS Github
BASE_URL = "https://github.com/ersilia-os/"
RAW_CONTENT_URL = "https://raw.githubusercontent.com/ersilia-os/{model}/main/"
REPO_API_URL = "https://api.github.com/repos/ersilia-os/{model}/contents"
USER_AGENT = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36"