"""Fetch model"""

from ... import ErsiliaBase


class ModelFetcher(ErsiliaBase):
    def __init__(self):
        pass

    def _get_repo(self):
        pass

    def _pack(self):
        pass

    def fetch(self, model_id):
        self.model_id = model_id
        if self.overwrite:
            pass
        ms = ModelStatus()
        if not ms.is_downloaded(self.model_id):
            self._get_repo()
            self._get_model_parameters()
            self._pack()
        if self.dockerize:
            if not ms.is_docker(self.model_id):
                self._dockerize()
        if self.pip:
            if not ms.is_pip(self.model_id):
                self._pip_install()
