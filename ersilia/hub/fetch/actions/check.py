from . import BaseAction
from ....serve.autoservice import AutoService


class ModelChecker(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )

    def check(self):
        self.logger.debug("Checking that autoservice works")
        AutoService(self.model_id)
