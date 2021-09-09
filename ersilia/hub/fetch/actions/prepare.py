from . import BaseAction
from ...bundle.status import ModelStatus
from ...delete.delete import ModelFullDeleter


class ModelPreparer(BaseAction):
    def __init__(self, model_id, overwrite, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.overwrite = overwrite
        self.status = ModelStatus(config_json=self.config_json)
        self.deleter = ModelFullDeleter(config_json=self.config_json)

    def prepare(self):
        if self.status.is_downloaded(self.model_id):
            if not self.overwrite:
                return
        self.deleter.delete(self.model_id)
