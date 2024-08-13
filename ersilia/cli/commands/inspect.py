from ...publish.inspect import ModelInspector
import click
from . import ersilia_cli
import json


def inspect_cmd():

    @ersilia_cli.command(short_help="Inspect model", help="Inspect model structure and metadata")
    @click.argument("model", type=click.STRING)
    def inspect(model):

        inspector = ModelInspector(model)
        value = {
             'is_github_url_available': inspector.checkRepoExists(0),
             'is_github_url_available_details': inspector.checkRepoExists(1),
             'metadata_complete': inspector.metadataComplete(0),
             'metadata_complete_details': inspector.metadataComplete(1),
             'folder_structure_complete': inspector.folderStructureComplete(0),
             'folder_structure_complete_details': inspector.folderStructureComplete(1),
             'docker_check': inspector.validateDependicies(0),
             'docker_check_details': inspector.validateDependicies(1),
             'computational_performance_tracking': inspector.computationalPerformance(0),
             'computational_performance_tracking_details': inspector.computationalPerformance(1),
             'extra_files_check': inspector.noExcessFiles(0),
             'extra_files_check_details': inspector.noExcessFiles(1),
        }
        
        print(json.dumps(value))
        return json.dumps(value)