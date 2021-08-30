from ..serve.autoservice import AutoService


class ErsiliaModel(AutoService):
    def __init__(self, model_id, config_json=None, overwrite=False):
        self.overwrite = overwrite
        self.config_json = config_json
        self.model_id = model_id
        AutoService.__init__(self, model_id, config_json=config_json)
