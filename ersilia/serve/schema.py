import os
import json
import numpy as np

from ..default import API_SCHEMA_FILE
from .. import ErsiliaBase


class ApiSchema(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.schema_file = os.path.join(
            self._model_path(self.model_id), API_SCHEMA_FILE
        )
        if not os.path.exists(self.schema_file):
            self.logger.debug("Schema not yet available")
        else:
            self.logger.debug("Schema available in {0}".format(self.schema_file))
        self._array_types = set(
            ["array", "mixed_array", "string_array", "numeric_array"]
        )
        self._h5_serializable_types = set(["numeric", "numeric_array", "array"])

    def _features(self, o):
        if o["meta"] is not None:
            return o["meta"]
        if o["type"] in self._array_types:
            shape = o["shape"]
        else:
            return None
        if len(shape) == 1:  # array
            n = shape[0]
            chars = len(str(n))
            names = []
            for i in range(n):
                i = str(i).zfill(chars)
                names += ["f{0}".format(i)]
            return names
        if len(shape) == 2:  # image
            n = shape[0]
            m = shape[1]
            names = []
            for i in range(n):
                names_ = []
                for j in range(m):
                    names_ += ["f{0}-{1}".format(i, j)]
                names += [names_]
            return names
        if len(shape) == 3:  # image with channels
            n = shape[0]
            m = shape[1]
            l = shape[2]
            names = []
            for i in range(n):
                names_ = []
                for j in range(m):
                    names__ = []
                    for k in range(l):
                        names__ += ["f{0}-{1}-{2}".format(i, j, k)]
                    names_ += [names__]
                names += [names_]
            return names
        # TODO: Make this generalizable
        return None

    def isfile(self):
        return os.path.isfile(self.schema_file)

    def get(self):
        with open(self.schema_file) as f:
            data = json.load(f)
        for api, sc in data.items():
            for k, o in sc["output"].items():
                data[api]["output"][k]["meta"] = self._features(o)
        return data

    @property
    def schema(self):
        return self.get()

    def get_schema_by_api(self, api_name):
        return self.schema[api_name]

    def get_output_by_api(self, api_name):
        return self.schema[api_name]["output"]

    def is_h5_serializable(self, api_name):
        schema = self.get_output_by_api(api_name)
        for _, v in schema.items():
            if v["type"] not in self._h5_serializable_types:
                return False
        return True

    def get_meta_by_api(self, api_name):
        sc = self.schema[api_name]["output"]
        meta = {}
        for k, v in sc.items():
            meta[k] = v["meta"]
        return meta

    def get_meta(self):
        sc = self.schema
        meta = {}
        for api, _ in sc.items():
            meta_ = self.get_meta_by_api(api)
            meta[api] = meta_
        return meta

    def get_apis(self):
        return sorted(self.schema.keys())

    def empty_by_field(self, field):
        if field["type"] in self._array_types:
            shape = tuple(field["shape"])
            return np.full(shape, None).tolist()
        return None

    def empty_input_by_api(self, api_name):
        sc = self.schema[api_name]["input"]
        d = {}
        for k, v in sc.items():
            d[k] = self.empty_by_field(v)
        return d

    def empty_output_by_api(self, api_name):
        sc = self.schema[api_name]["output"]
        d = {}
        for k, v in sc.items():
            d[k] = self.empty_by_field(v)
        return d

    def empty_by_api(self, api_name):
        return {
            "input": self.empty_input_by_api(api_name),
            "output": self.empty_output_by_api(api_name),
        }

    def empty(self):
        d = {}
        for api_name in self.get_apis():
            d[api_name] = self.empty_by_api(api_name)
        return d
