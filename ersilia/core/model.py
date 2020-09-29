import importlib
from .base import ErsiliaBase
from ..hub.fetch import ModelFetcher


class ErsiliaModel(ErsiliaBase):

    def __init__(self, model_id, config_json=None, pip=True):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.pip = pip
        self._fetch()
        self._import()

    def _fetch(self, overwrite=True, local=False):
        mf = ModelFetcher(config_json=self.config_json, overwrite=overwrite, local=local)
        mf.fetch(self.model_id, pip=self.pip)

    def _import(self):
        Model = importlib.import_module(self.model_id, package=None)
        self.mdl = Model.load()

    def predict(self, inp):
        return self.mdl.predict(inp)

