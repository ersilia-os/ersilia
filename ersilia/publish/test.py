import os

from ..hub.content.card import RepoMetadataFile
from .. import ErsiliaBase

from ..default import METADATA_JSON_FILE


class LocalModelTester(ErsiliaBase):
    def __init__(self, model_id, repo_path, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.repo_path = os.path.abspath(repo_path)

    def check_metadata(self):
        json_file = os.path.join(self.repo_path, METADATA_JSON_FILE)
        rm = RepoMetadataFile()
        rm.read_information(json_path=json_file)
        self.logger.success("Metadata check successful")

    def check_fetch(self):
        pass  # TODO
        self.logger.success("Fetch check successful")


class RemoteModelTester(ErsiliaBase):
    pass  # TODO
