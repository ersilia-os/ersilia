from . import BaseAction
from ...bundle.status import ModelStatus
from ...delete.delete import ModelFullDeleter

from ....utils.exceptions_utils.delete_exceptions import ModelDeleteError
from .... import throw_ersilia_exception


class ModelPreparer(BaseAction):
    def __init__(self, model_id, overwrite, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.overwrite = overwrite
        self.status = ModelStatus(config_json=self.config_json)
        self.deleter = ModelFullDeleter(config_json=self.config_json)

    @throw_ersilia_exception
    def prepare(self):
        if self.status.is_downloaded(self.model_id):
            if not self.overwrite:
                return
        try:
            self.deleter.delete(self.model_id)
        except:
            raise ModelDeleteError(model=self.model_id)
