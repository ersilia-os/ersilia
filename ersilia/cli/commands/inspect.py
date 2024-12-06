import click
import json
from . import ersilia_cli
from ...publish.inspect import ModelInspector
from ...cli import echo
from ... import ErsiliaBase

# Dont put an echo here: the json output will be used in test cmd
# try: logger
class InspectCommand(ErsiliaBase):
    def __init__(
        self, 
        config_json=None, 
        credentials_json=None
    ):
        super().__init__(config_json, credentials_json)

    def _perform_checks(self, inspector, dir):
      
        checks = {
            "repo_result": inspector.check_repo_exists if dir is None else None,
            "metadata_result": inspector.check_complete_metadata if dir is None else None,
            "folder_structure_result": inspector.check_complete_folder_structure,
            "dependencies_result": inspector.check_dependencies_are_valid,
            "performance_result": inspector.check_computational_performance,
            "extra_files_result": inspector.check_no_extra_files if dir is None else None,
        }

        return {key: check() for key, check in checks.items() if check}

    def _format_output(self, results, dir):
     
        output = {}
        mapping = {
            "repo_result": (
                "is_github_url_available", 
                "is_github_url_available_details"
            ),
            "metadata_result": (
                "complete_metadata", 
                "complete_metadata_details"
            ),
            "folder_structure_result": (
                "complete_folder_structure", 
                "complete_folder_structure_details"
            ),
            "dependencies_result": (
                "docker_check", 
                "docker_check_details"
            ),
            "performance_result": (
                "computational_performance_tracking", 
                "computational_performance_tracking_details"
            ),
            "extra_files_result": (
                "extra_files_check", 
                "extra_files_check_details"
            ),
        }

        for key, result in results.items():
            success_key, details_key = mapping[key]
            output[success_key] = result.success
            output[details_key] = result.details

        return output

    def inspect(self, model, dir):
        inspector = ModelInspector(model, dir)
        results = self._perform_checks(inspector, dir)
        output = self._format_output(results, dir)
        echo(json.dumps(output, indent=2))
        return json.dumps(output)


def inspect_cmd():
    @ersilia_cli.command(
        short_help="Inspect model",
        help="Inspect model structure and metadata."
    )
    @click.argument("model", type=click.STRING)
    @click.option(
        "-d",
        "--dir",
        "dir",
        required=False,
        default=None,
        type=click.STRING,
    )
    def inspect(model, dir):
        cmd = InspectCommand()
        return cmd.inspect(model, dir)
