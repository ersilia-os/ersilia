import json
import boto3

from ...default import ERSILIA_MODEL_HUB_S3_BUCKET, MODELS_JSON
from ... import ErsiliaBase

import requests


class JsonModelsInterface(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.json_file_name = MODELS_JSON
        self.url = f"https://{ERSILIA_MODEL_HUB_S3_BUCKET}.s3.eu-central-1.amazonaws.com/{MODELS_JSON}"

    def _read_json_file(self):
        response = requests.get(self.url)
        models_list = response.json()
        return models_list

    def items(self):
        models = self._read_json_file()
        for mdl in models:
            yield mdl

    def items_all(self):
        models = self._read_json_file()
        return models
