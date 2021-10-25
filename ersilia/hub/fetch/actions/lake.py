from ....utils.dvc import DVCFetcher

from . import BaseAction


class LakeGetter(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.model_dest = self._model_path(model_id)
        self.dvc_fetcher = DVCFetcher(self.model_dest)

    def get(self):
        self.dvc_fetcher.get_data()
