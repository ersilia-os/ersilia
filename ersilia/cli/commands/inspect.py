from ...publish.inspect import ModelInspector
import click
from . import ersilia_cli
import json
from ... import ErsiliaBase


class InspectCommand(ErsiliaBase):
    def __init__(self, config_json=None, credentials_json=None):
        super().__init__(config_json, credentials_json)

    def inspect(self, model):
        inspector = ModelInspector(model)

        check_repo_exists_result = inspector.check_repo_exists()
        check_complete_metadata_result = inspector.check_complete_metadata()
        check_complete_folder_structure_result = (
            inspector.check_complete_folder_structure()
        )
        check_dependencies_are_valid_result = inspector.check_dependencies_are_valid()
        check_comptuational_performance_result = (
            inspector.check_comptuational_performance()
        )
        check_no_extra_files_result = inspector.check_no_extra_files()

        value = {
            "is_github_url_available": check_repo_exists_result.success,
            "is_github_url_available_details": check_repo_exists_result.details,
            "complete_metadata": check_complete_metadata_result.success,
            "complete_metadata_details": check_complete_metadata_result.details,
            "complete_folder_structure": check_complete_folder_structure_result.success,
            "complete_folder_structure_details": check_complete_folder_structure_result.details,
            "docker_check": check_dependencies_are_valid_result.success,
            "docker_check_details": check_dependencies_are_valid_result.details,
            "computational_performance_tracking": check_comptuational_performance_result.success,
            "computational_performance_tracking_details": check_comptuational_performance_result.details,
            "extra_files_check": check_no_extra_files_result.success,
            "extra_files_check_details": check_no_extra_files_result.details,
        }

        self.logger.debug(json.dumps(value, indent=2))
        return json.dumps(value)


def inspect_cmd():
    @ersilia_cli.command(
        short_help="Inspect model", help="Inspect model structure and metadata"
    )
    @click.argument("model", type=click.STRING)
    def inspect(model):
        cmd = InspectCommand()
        return cmd.inspect(model)
