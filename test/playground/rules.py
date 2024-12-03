import json
import csv
from pathlib import Path

RULE_REGISTRY = {}


class CommandRule:
    def check(self, *args, **kwargs):
        raise NotImplementedError(
            "Each rule must implement a check method."
        )


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
        actual_status = Path(folder_path).exists() and any(
            Path(folder_path).iterdir()
        )
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
            actual_status = (
                f'"docker_hub": {str(expected_status).lower()}'
                in content
            )
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
    def __init__(self):
        pass

    def check(self, file_path, expected_status):
        if not Path(file_path).exists():
            raise FileNotFoundError(
                f"File {file_path} does not exist."
            )

        file_extension = Path(file_path).suffix.lower()
        if file_extension not in [".json", ".csv"]:
            raise ValueError(
                f"Unsupported file type: {file_extension}. Only JSON and CSV are supported."
            )

        if file_extension == ".json":
            actual_status = self._check_json_content(file_path)
        elif file_extension == ".csv":
            actual_status = self._check_csv_content(file_path)
        else:
            raise ValueError(
                f"Unexpected error occurred with file extension: {file_extension}"
            )

        if actual_status != expected_status:
            raise AssertionError(
                f"Expectation failed for FileContentCheckRule: "
                f"Expected file content to be '{expected_status}', "
                f"but it was '{actual_status}'."
            )

        return {
            "name": f"File content check at {file_path}",
            "status": actual_status,
        }

    def _check_json_content(self, file_path):
        """Checks the content of a JSON file."""
        with open(file_path, "r") as f:
            try:
                content = json.load(f)
                return "not null" if content else "null"
            except json.JSONDecodeError as e:
                raise ValueError(
                    f"Invalid JSON content in file {file_path}: {e}"
                )

    def _check_csv_content(self, file_path):
        """Checks the content of a CSV file."""
        with open(file_path, "r") as f:
            reader = csv.reader(f)
            try:
                rows = list(reader)
                return "not null" if len(rows) > 1 else "null"
            except csv.Error as e:
                raise ValueError(
                    f"Invalid CSV content in file {file_path}: {e}"
                )


def get_rule(rule_name, *args, **kwargs):
    rule_class = RULE_REGISTRY.get(rule_name)
    if not rule_class:
        raise ValueError(f"Rule '{rule_name}' is not registered.")
    return rule_class().check(*args, **kwargs)
