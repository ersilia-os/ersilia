import importlib
from .base import ErsiliaBase
from ..hub.fetch import ModelFetcher


class ErsiliaModel(ErsiliaBase):

    def __init__(self, model_id, config_json=None, overwrite=True, pip=True, local=False):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.overwrite = overwrite
        self.local = local
        self.model_id = model_id
        self.pip = pip
        self._fetch()
        self._import()

    def _fetch(self):
        mf = ModelFetcher(config_json=self.config_json, overwrite=self.overwrite, local=self.local)
        mf.fetch(self.model_id, pip=self.pip)

    def _import(self):
        Model = importlib.import_module(self.model_id, package=None)
        self.mdl = Model.load()

    def predict(self, inp):
        return self.mdl.predict(inp)
