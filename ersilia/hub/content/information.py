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


class Information(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.repository_folder = os.path.join(
            self._get_bento_location(model_id=self.model_id)
        )
        self.dest_folder = os.path.join(self._model_path(model_id=model_id))

    def _get_pack_mode(self):
        pack_mode_file = os.path.join(self.dest_folder, PACKMODE_FILE)
        if os.path.exists(pack_mode_file):
            with open(pack_mode_file, "r") as f:
                return f.read().rstrip()
        else:
            return None

    def _get_service_class(self):
        service_class_file = os.path.join(self.repository_folder, SERVICE_CLASS_FILE)
        if os.path.exists(service_class_file):
            with open(service_class_file, "r") as f:
                return f.read().rstrip()
        else:
            return None

    def _get_api_schema(self):
        api_schema_file = os.path.join(self.dest_folder, API_SCHEMA_FILE)
        if os.path.exists(api_schema_file):
            with open(api_schema_file, "r") as f:
                return json.load(f)
        else:
            return None

    def _get_size(self):
        size_file = os.path.join(self.dest_folder, MODEL_SIZE_FILE)
        if os.path.exists(size_file):
            with open(size_file, "r") as f:
                return json.load(f)
        else:
            return None

    def _get_metadata(self):
        metadata_file = os.path.join(self.dest_folder, METADATA_JSON_FILE)
        if os.path.exists(metadata_file):
            with open(metadata_file, "r") as f:
                return json.load(f)
        else:
            return None

    def _get_card(self):
        card_file = os.path.join(self.dest_folder, CARD_FILE)
        if os.path.exists(card_file):
            with open(card_file, "r") as f:
                return json.load(f)
        else:
            return None

    def _get_apis_list(self):
        apis_list_file = os.path.join(self.dest_folder, CARD_FILE)
        if os.path.exists(apis_list_file):
            with open(os.path.join(self.repository_folder, APIS_LIST_FILE), "r") as f:
                return [x.rstrip() for x in f.readlines()]
        else:
            return None

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
