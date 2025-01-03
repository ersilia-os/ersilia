import json
import os

import numpy as np

from .. import ErsiliaBase
from ..default import API_SCHEMA_FILE


class ApiSchema(ErsiliaBase):
    """
    Class to manage and interact with the API schema for a given model.

    An API schema defines the structure of the input and output data for the APIs of a model.
    It includes information about the types and shapes of the data fields, as well as any metadata
    associated with them. This class provides methods to retrieve and manipulate the schema, check
    for HDF5 serializability, and generate empty representations of the input and output data.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    config_json : dict
        Configuration in JSON format.
    """

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
            if n is None:
                return None
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

    def isfile(self) -> bool:
        """
        Check if the schema file exists.

        Returns
        -------
        bool
            True if the schema file exists, False otherwise.
        """
        return os.path.isfile(self.schema_file)

    def get(self) -> dict:
        """
        Get the schema data.

        Returns
        -------
        dict
            The schema data.
        """
        with open(self.schema_file) as f:
            data = json.load(f)
        for api, sc in data.items():
            for k, o in sc["output"].items():
                data[api]["output"][k]["meta"] = self._features(o)
        return data

    @property
    def schema(self) -> dict:
        """
        Get the schema data.

        Returns
        -------
        dict
            The schema data.
        """
        return self.get()

    def get_schema_by_api(self, api_name: str) -> dict:
        """
        Get the schema for a specific API.

        Parameters
        ----------
        api_name : str
            The name of the API.

        Returns
        -------
        dict
            The schema for the specified API.
        """
        return self.schema[api_name]

    def get_output_by_api(self, api_name: str) -> dict:
        """
        Get the output schema for a specific API.

        Parameters
        ----------
        api_name : str
            The name of the API.

        Returns
        -------
        dict
            The output schema for the specified API.
        """
        return self.schema[api_name]["output"]

    def is_h5_serializable(self, api_name: str) -> bool:
        """
        Check if the output schema for a specific API is HDF5 serializable.

        Parameters
        ----------
        api_name : str
            The name of the API.

        Returns
        -------
        bool
            True if the output schema is HDF5 serializable, False otherwise.
        """
        schema = self.get_output_by_api(api_name)
        for _, v in schema.items():
            if v["type"] not in self._h5_serializable_types:
                return False
        return True

    def get_meta_by_api(self, api_name: str) -> dict:
        """
        Get the metadata for a specific API.

        Parameters
        ----------
        api_name : str
            The name of the API.

        Returns
        -------
        dict
            The metadata for the specified API.
        """
        sc = self.schema[api_name]["output"]
        meta = {}
        for k, v in sc.items():
            meta[k] = v["meta"]
        return meta

    def get_meta(self) -> dict:
        """
        Get the metadata for all APIs.

        Returns
        -------
        dict
            The metadata for all APIs.
        """
        sc = self.schema
        meta = {}
        for api, _ in sc.items():
            meta_ = self.get_meta_by_api(api)
            meta[api] = meta_
        return meta

    def get_apis(self) -> list:
        """
        Get the list of APIs available for the model.

        Returns
        -------
        list
            List of API names.
        """
        return sorted(self.schema.keys())

    def empty_by_field(self, field: dict) -> list:
        """
        Get an empty representation for a specific field.

        Parameters
        ----------
        field : dict
            The field schema.

        Returns
        -------
        list
            An empty representation of the field.
        """
        if field["type"] in self._array_types:
            shape = tuple(field["shape"])
            return np.full(shape, None).tolist()
        return None

    def empty_input_by_api(self, api_name: str) -> dict:
        """
        Get an empty input representation for a specific API.

        Parameters
        ----------
        api_name : str
            The name of the API.

        Returns
        -------
        dict
            An empty input representation for the specified API.
        """
        sc = self.schema[api_name]["input"]
        d = {}
        for k, v in sc.items():
            d[k] = self.empty_by_field(v)
        return d

    def empty_output_by_api(self, api_name: str) -> dict:
        """
        Get an empty output representation for a specific API.

        Parameters
        ----------
        api_name : str
            The name of the API.

        Returns
        -------
        dict
            An empty output representation for the specified API.
        """
        sc = self.schema[api_name]["output"]
        d = {}
        for k, v in sc.items():
            d[k] = self.empty_by_field(v)
        return d

    def empty_by_api(self, api_name: str) -> dict:
        """
        Get an empty input and output representation for a specific API.

        Parameters
        ----------
        api_name : str
            The name of the API.

        Returns
        -------
        dict
            An empty input and output representation for the specified API.
        """
        return {
            "input": self.empty_input_by_api(api_name),
            "output": self.empty_output_by_api(api_name),
        }

    def empty(self) -> dict:
        """
        Get an empty input and output representation for all APIs.

        Returns
        -------
        dict
            An empty input and output representation for all APIs.
        """
        d = {}
        for api_name in self.get_apis():
            d[api_name] = self.empty_by_api(api_name)
        return d
