import os
import json

from ... import ErsiliaBase
from ...default import (
    PACKMODE_FILE,
    API_SCHEMA_FILE,
    MODEL_SIZE_FILE,
    METADATA_JSON_FILE,
    CARD_FILE,
    SERVICE_CLASS_FILE,
    APIS_LIST_FILE,
)
from ..fetch import STATUS_FILE


class Information(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.repository_folder = os.path.join(
            self._get_bento_location(model_id=self.model_id)
        )
        self.dest_folder = os.path.join(self._model_path(model_id=model_id))

    def _get_pack_mode(self):
        with open(os.path.join(self.dest_folder, PACKMODE_FILE), "r") as f:
            return f.read().rstrip()

    def _get_service_class(self):
        with open(os.path.join(self.repository_folder, SERVICE_CLASS_FILE), "r") as f:
            return f.read().rstrip()

    def _get_api_schema(self):
        with open(os.path.join(self.dest_folder, API_SCHEMA_FILE), "r") as f:
            return json.load(f)

    def _get_size(self):
        with open(os.path.join(self.dest_folder, MODEL_SIZE_FILE), "r") as f:
            return json.load(f)

    def _get_metadata(self):
        with open(os.path.join(self.dest_folder, METADATA_JSON_FILE), "r") as f:
            return json.load(f)

    def _get_card(self):
        with open(os.path.join(self.dest_folder, CARD_FILE), "r") as f:
            return json.load(f)

    def _get_apis_list(self):
        with open(os.path.join(self.repository_folder, APIS_LIST_FILE), "r") as f:
            return [x.rstrip() for x in f.readlines()]

    def get(self):
        data = {
            "pack_mode": self._get_pack_mode(),
            "service_class": self._get_service_class(),
            "apis_list": self._get_apis_list(),
            "api_schema": self._get_api_schema(),
            "size": self._get_size(),
            "metadata": self._get_metadata(),
            "card": self._get_card(),
        }
        return data
