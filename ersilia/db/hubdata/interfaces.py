import requests

from ... import ErsiliaBase
from ...default import ERSILIA_MODEL_HUB_S3_BUCKET, MODELS_JSON


class JsonModelsInterface(ErsiliaBase):
    """
    Interface for fetching model metadata from models stored in S3.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for initializing the interface.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.json_file_name = MODELS_JSON
        self.url = f"https://{ERSILIA_MODEL_HUB_S3_BUCKET}.s3.eu-central-1.amazonaws.com/{MODELS_JSON}"

    def _read_json_file(self):
        response = requests.get(self.url)
        models_list = response.json()
        return models_list

    def items(self):
        """
        Yields models from the JSON file one by one.

        Yields
        ------
        dict
            A model from the JSON file.
        """
        models = self._read_json_file()
        for mdl in models:
            yield mdl

    def items_all(self):
        """
        Retrieves all models from the JSON file.

        Returns
        -------
        list
            List of all models.
        """
        models = self._read_json_file()
        return models

    def identifier_exists(self, model_id):
        """
        Checks if a model identifier exists in the JSON file.

        Parameters
        ----------
        model_id : str
            Identifier of the model to check.

        Returns
        -------
        bool
            True if the identifier exists, False otherwise.
        """
        data = self._read_json_file()
        return any(item["Identifier"] == model_id for item in data)
