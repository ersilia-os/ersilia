from ....utils.dvc import DVCFetcher
from . import BaseAction


class LakeGetter(BaseAction):
    """
    Class to fetch data from precalculated data from the DVC repository.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    config_json : dict
        Configuration settings for the model.
    """

    def __init__(self, model_id: str, config_json: dict):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.model_dest = self._model_path(model_id)
        self.dvc_fetcher = DVCFetcher(self.model_dest)

    def get(self):
        """
        Fetch data using DVCFetcher.
        """
        self.dvc_fetcher.get_data()
