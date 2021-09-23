import os
import json
import numpy as np

from ..default import API_SCHEMA_FILE
from .. import ErsiliaBase


class ApiSchema(ErsiliaBase):

    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.schema = self.get()

    def get(self):
        with open(
            os.path.join(self._model_path(self.model_id), API_SCHEMA_FILE), "r"
        ) as f:
            return json.load(f)

    def get_apis(self):
        return sorted(self.schema.keys())

    def empty_by_field(self, field):
        if field["type"] == "array":
            shape = tuple(field["shape"])
            return np.full(shape, None).tolist()
        return None

    def empty_input_by_api(self, api_name):
        sc = self.schema[api_name]["input"]
        d = {}
        for k,v in sc.items():
            d[k] = self.empty_by_field(v)
        return d

    def empty_output_by_api(self, api_name):
        sc = self.schema[api_name]["output"]
        d = {}
        for k,v in sc.items():
            d[k] = self.empty_by_field(v)
        return d

    def empty_by_api(self, api_name):
        return {"input": self.empty_input_by_api(api_name),
                "output": self.empty_output_by_api(api_name)}

    def empty(self):
        d = {}
        for api_name in self.get_apis():
            d[api_name] = self.empty_by_api(api_name)
        return d
