# Description: This file contains the class that fetches data from the data lake.

from . import BaseAction


class LakeGetter(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.model_dest = self._model_path(model_id)
        # TODO Initialize connection to the data lake that we will use

    def get(self):
        # TODO Fetch the data from the data lake
        raise NotImplementedError("This feature is not yet implemented.")

