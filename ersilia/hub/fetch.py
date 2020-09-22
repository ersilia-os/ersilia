"""Fetch model and pack it to a BentoML bundle"""

import os
import runpy
from ersilia import Config
from ersilia.utils.download import GitHubDownloader


class ModelFetcher(object):

    def __init__(self):
        self.data = Config.URL.DATA
        self.dest = Config.URL.DEST
        self.token = Config.HUB.TOKEN
        self.org = Config.HUB.ORG
        self.repo = Config.HUB.REPO
        self.tag = Config.HUB.TAG
        self.down = GitHubDownloader(self.token)

    def _model_path(self, model_id):
        folder = os.path.join(self.dest, "models", model_id)
        return folder

    def get_repo(self, model_id):
        folder = self._model_path(model_id)
        if os.path.exists(folder):
            pass
        else:
            repo_folder = "models/"+model_id
            self.down.fetch(repo_folder, org=self.org, repo=self.repo, tag=self.tag)

    def get_data(self, model_id):
        folder = self._model_path(model_id)


    def pack(self, model_id):
        folder = self._model_path(model_id)
        pack_script = os.path.join(folder, Config.HUB.PACK_SCRIPT)
        runpy.run_path(path_name=pack_script)

    def fetch(self, model_id):
        self.get_repo(model_id)
        self.get_data(model_id)
        self.pack(model_id)