from . import BaseAction
from .... import ErsiliaModel


class ModelSniffer(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.model = ErsiliaModel(model_id, config_json=config_json, overwrite=False)

    def _size(self):
        pass

    def _speed(self):
        pass

    def _output_schema(self):
        pass

    def sniff(self):
        self.logger.debug("Sniffing model")
        AutoService(self.model_id)
