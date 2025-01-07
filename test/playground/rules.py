import csv
import json
import math
import yaml
from ersilia.utils.hdf5 import Hdf5DataLoader
from pathlib import Path

RULE_REGISTRY = {}

config_path = Path("config.yml")
config = yaml.safe_load(config_path.read_text())


class CommandRule:
    def check(self, *args, **kwargs):
        raise NotImplementedError("Each rule must implement a check method.")


def register_rule(name):
    def decorator(cls):
        RULE_REGISTRY[name] = cls
        return cls

    return decorator


@register_rule("folder_exists")
class FolderExistsRule(CommandRule):
    def __init__(self):
        pass

    def check(self, folder_path, expected_status):
        actual_status = Path(folder_path).exists() and any(Path(folder_path).iterdir())
        if actual_status != expected_status:
            raise AssertionError(
                f"Expectation failed for FolderExistsRule: "
                f"Expected folder to {'exist' if expected_status else 'not exist'}, "
                f"but it {'exists' if actual_status else 'does not exist'}."
            )
        return {
            "name": f"Folder exists at {folder_path}",
            "status": actual_status,
        }


@register_rule("file_exists")
class FileExistsRule(CommandRule):
    def __init__(self):
        pass

    def check(self, file_path, expected_status):
        actual_status = Path(file_path).exists()
        if actual_status != expected_status:
            raise AssertionError(
                f"Expectation failed for FileExistsRule: "
                f"Expected file to {'exist' if expected_status else 'not exist'}, "
                f"but it {'exists' if actual_status else 'does not exist'}."
            )
        return {
            "name": f"File exists at {file_path}",
            "status": actual_status,
        }


@register_rule("dockerhub_status")
class DockerHubStatusRule(CommandRule):
    def __init__(self):
        pass

    def check(self, expected_status, dest_path):
        dockerhub_file = Path(dest_path) / "from_dockerhub.json"
        if dockerhub_file.exists():
            with open(dockerhub_file, "r") as f:
                content = f.read()
            actual_status = f'"docker_hub": {str(expected_status).lower()}' in content
        else:
            actual_status = False

        if actual_status != expected_status:
            raise AssertionError(
                f"Expectation failed for DockerHubStatusRule: "
                f"Expected DockerHub status to be {expected_status}, but it was {actual_status}."
            )
        return {
            "name": f"DockerHub status is {actual_status}",
            "status": actual_status,
        }


@register_rule("file_content_check")
class FileContentCheckRule(CommandRule):
    # Supported data structures
    ds = {
        "Single": lambda x: isinstance(x, list) and len(x) == 1,
        "List": lambda x: isinstance(x, list)
        and len(x) > 1
        and all(isinstance(item, (int, float)) for item in x),
        "Flexible List": lambda x: isinstance(x, list)
        and all(isinstance(item, (str, int, float)) for item in x),
        "Matrix": lambda x: isinstance(x, list)
        and all(
            isinstance(row, list)
            and all(isinstance(item, (int, float)) for item in row)
            for row in x
        ),
        "Serializable Object": lambda x: isinstance(x, dict),
    }

    def __init__(self):
        pass

    def check(self, file_path, expected_input):
        if not Path(file_path).exists():
            raise FileNotFoundError(f"File {file_path} does not exist.")

        file_extension = Path(file_path).suffix.lower()
        if file_extension not in [".json", ".csv", ".h5"]:
            raise ValueError(
                f"Unsupported file type: {file_extension}. Only JSON, CSV, and HDF5 are supported."
            )

        if file_extension == ".json":
            actual_status, file_length = self._check_json_content(file_path)
        elif file_extension == ".csv":
            actual_status, file_length = self._check_csv_content(file_path)
        elif file_extension == ".h5":
            actual_status, file_length = self._check_h5_content(file_path)

        expected_length = self._determine_expected_length(expected_input)

        if file_length != expected_length:
            raise AssertionError(
                f"Expectation failed for FileContentCheckRule: "
                f"Expected length '{expected_length}', but found '{file_length}' in the file."
            )

        return {
            "name": f"File content check at {file_path}",
            "status": actual_status,
            "length": file_length,
        }

    def _check_json_content(self, file_path):
        with open(file_path, "r") as f:
            try:
                content = json.load(f)
                for row in content:
                    row = row.get("output", {}).get("outcome")
                    self._validate_json(row)
                length = len(content) if isinstance(content, list) else 1
                return "not null", length
            except json.JSONDecodeError as e:
                raise ValueError(f"Invalid JSON content in file {file_path}: {e}")

    def _check_csv_content(self, file_path):
        with open(file_path, "r") as f:
            reader = csv.reader(f)
            rows = list(reader)[1:]
            for index, row in enumerate(rows):
                if any(self._is_invalid_value(cell) for cell in row):
                    raise ValueError(
                        f"Invalid value found in row {index + 1} of CSV file."
                    )
            return "not null", len(rows)

    def _check_h5_content(self, file_path):
        try:
            loader = Hdf5DataLoader()
            loader.load(file_path)
            content = loader.values or loader.keys or loader.inputs or loader.features
            if not content:
                raise ValueError(f"Empty content in HDF5 file {file_path}.")
            content = (
                content[1:]
                if isinstance(content, list) and len(content) > 1
                else content
            )
            for index, row in enumerate(content):
                if any(self._is_invalid_value(cell) for cell in row):
                    raise ValueError(
                        f"Invalid value found in row {index + 1} of HDF5 file."
                    )
            return "not null", len(content)
        except (OSError, KeyError, ValueError) as e:
            raise ValueError(f"Invalid HDF5 content in file {file_path}: {e}")

    def _determine_expected_length(self, expected_input):
        if "str" in expected_input:
            return 1
        else:
            return config["number_of_input_samples"]

    def _validate_json(self, outcome):
        for shape, validator in self.ds.items():
            if validator(outcome):
                print(f"Outcome type identified as: {shape}")

                if shape in ["Single"]:
                    if self._is_invalid_value(outcome[0]):
                        raise ValueError(
                            f"Invalid value '{outcome[0]}' in outcome of type '{shape}'."
                        )
                elif shape in ["List", "Flexible List"]:
                    for item in outcome:
                        if self._is_invalid_value(item):
                            raise ValueError(
                                f"Invalid value '{item}' in outcome of type '{shape}'."
                            )
                elif shape == "Matrix":
                    for row in outcome:
                        for item in row:
                            if self._is_invalid_value(item):
                                raise ValueError(
                                    f"Invalid value '{item}' in outcome of type '{shape}'."
                                )
                elif shape == "Serializable Object":
                    for key, value in outcome.items():
                        if self._is_invalid_value(value):
                            raise ValueError(
                                f"Invalid value '{value}' found for key '{key}' in outcome of type '{shape}'."
                            )
                return

        raise TypeError("Unknown outcome structure.")

    def _is_invalid_value(self, item):
        if item in [None, "", "null"]:
            return True
        if isinstance(item, float) and math.isnan(item):
            return True
        return False


def get_rule(rule_name, *args, **kwargs):
    rule_class = RULE_REGISTRY.get(rule_name)
    if not rule_class:
        raise ValueError(f"Rule '{rule_name}' is not registered.")
    return rule_class().check(*args, **kwargs)
