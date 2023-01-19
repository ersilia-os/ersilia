import os
import tempfile

from ..hub.content.card import ReadmeMetadata, AirtableMetadata, RepoMetadataFile
from ..utils.terminal import run_command
from .. import ErsiliaBase

from ..default import GITHUB_ORG


class ReadmeUpdater(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        self.model_id = model_id
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-os")
        self.cwd = os.getcwd()
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)

    def _git_clone(self):
        run_command(
            "cd {0}; git clone https://github.com/{1}/{2}; cd {3}".format(
                self.tmp_folder, GITHUB_ORG, self.model_id, self.cwd
            )
        )

    def _git_push(self):
        run_command(
            'cd {0}/{1}; git add .; git commit -m "Updating README file from AirTable [skip ci]"; git push; cd {2}'.format(
                self.tmp_folder, self.model_id, self.cwd
            )
        )

    def update(self):
        self._git_clone()
        rm = ReadmeMetadata(model_id=self.model_id)
        bi = rm.read_information()
        tmp_file = os.path.join(self.tmp_folder, self.model_id, "README.md")
        rm.write_information(data=bi, readme_path=tmp_file)
        self._git_push()


class JsonUpdater(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        self.model_id = model_id
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-os")
        self.cwd = os.getcwd()
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)

    def _git_clone(self):
        run_command(
            "cd {0}; git clone https://github.com/{1}/{2}; cd {3}".format(
                self.tmp_folder, GITHUB_ORG, self.model_id, self.cwd
            )
        )

    def _git_push(self):
        run_command(
            'cd {0}/{1}; git add .; git commit -m "Updating metadata file from AirTable [skip ci]"; git push; cd {2}'.format(
                self.tmp_folder, self.model_id, self.cwd
            )
        )

    def update(self):
        self._git_clone()
        ai = AirtableMetadata(model_id=self.model_id)
        data = ai.read_information()
        tmp_file = os.path.join(self.tmp_folder, self.model_id, "metadata.json")
        rm = RepoMetadataFile(model_id=self.model_id)
        rm.write_information(data, tmp_file)
        self._git_push()
