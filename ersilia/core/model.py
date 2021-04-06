import importlib
from .base import ErsiliaBase
from ..hub.fetch import ModelFetcher
from ..serve.autoservice import AutoService


class ErsiliaModel(AutoService):

    def __init__(self, model_id, config_json=None, overwrite=False):
        self.overwrite = overwrite
        self.config_json = config_json
        self.model_id = model_id
        self._fetch()
        AutoService.__init__(self, model_id, config_json=config_json)

    def _fetch(self):
        mf = ModelFetcher(config_json=self.config_json, overwrite=self.overwrite)
        mf.fetch(self.model_id)
