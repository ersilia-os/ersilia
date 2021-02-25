import importlib
from .base import ErsiliaBase
from ..hub.fetch import ModelFetcher


class ErsiliaModel(ErsiliaBase):

    def __init__(self, model_id, config_json=None, overwrite=False, pip=True, local=True):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.overwrite = overwrite
        self.local = local
        self.model_id = model_id
        self.pip = pip
        self._fetch()
        self._pip_model = None

    def _fetch(self):
        mf = ModelFetcher(config_json=self.config_json, overwrite=self.overwrite, local=self.local)
        mf.fetch(self.model_id, pip=self.pip)

    def _pip_import(self):
        try:
            Model = importlib.import_module(self.model_id, package=None)
            return Model.load()
        except ModuleNotFoundError:
            return None

    def predict_with_pip(self, inp):
        if not self._pip_model:
            self._pip_model = self._pip_import()
        if not self._pip_model:
            pass

    def predict(self, inp):
        if type(inp) == str:
            is_single = True
        else:
            is_single = False
        if is_single:
            inp = [inp]
        pred = self.mdl.predict(inp)
        if is_single:
            return pred[0]
        else:
            return pred

    def serve(self):
        pass

    def close(self):
        pass
