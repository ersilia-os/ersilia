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
        value = {
            'is_github_url_available': inspector.check_repo_exists(0),
            'is_github_url_available_details': inspector.check_repo_exists(1),
            'metadata_complete': inspector.metadata_complete(0),
            'metadata_complete_details': inspector.metadata_complete(1),
            'folder_structure_complete': inspector.folder_structure_complete(0),
            'folder_structure_complete_details': inspector.folder_structure_complete(1),
            'docker_check': inspector.validate_dependencies(0),
            'docker_check_details': inspector.validate_dependencies(1),
            'computational_performance_tracking': inspector.computational_performance(0),
            'computational_performance_tracking_details': inspector.computational_performance(1),
            'extra_files_check': inspector.no_excess_files(0),
            'extra_files_check_details': inspector.no_excess_files(1),
        }
        
        self.logger.debug(json.dumps(value, indent=2))
        return json.dumps(value)

def inspect_cmd():
    @ersilia_cli.command(short_help="Inspect model", help="Inspect model structure and metadata")
    @click.argument("model", type=click.STRING)
    def inspect(model):
        cmd = InspectCommand()
        return cmd.inspect(model)