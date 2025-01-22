import json
import os

import requests

from .... import EOS, ErsiliaBase
from ....default import API_SCHEMA_FILE, INFORMATION_FILE, IS_FETCHED_FROM_HOSTED_FILE
from ....serve.services import HostedService
from ...fetch import ModelURLResolver
from .. import STATUS_FILE
from ..register.register import ModelRegisterer


class ModelHostedFetcher(ErsiliaBase):
    """
    A class used to fetch models from a hosted URL.

    Attributes
    ----------
    url : str
        URL where the model is hosted.
    config_json : dict
        Configuration settings in JSON format.

    Methods
    -------
    is_available(model_id)
        Check if the model is available at the hosted URL.
    write_apis(model_id)
        Write APIs for the model.
    get_information(model_id)
        Get information for the model.
    get_metadata(model_id)
        Get metadata for the model.
    write_status(model_id)
        Write the status file for the model.
    fetch(model_id)
        Fetch the model from the hosted URL.
    """

    def __init__(self, url: str, config_json: dict = None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.logger.debug("Initialized with URL: {0}".format(url))
        self.url = url

    def _is_available_known_url(self) -> bool:
        self.logger.debug("Checking if url {0} is reachable".format(self.url))
        try:
            response = requests.get(self.url, timeout=5)
            response.raise_for_status()
            return True
        except requests.exceptions.RequestException as err:
            self.logger.debug(
                "URL {0} is not reachable. Error: {1}".format(self.url, err)
            )
            return False

    def _is_available_unknown_url(self, model_id: str) -> bool:
        self.logger.debug(
            "Trying to find an available URL where the model is hosted using Models JSON"
        )
        mdl_url_resolver = ModelURLResolver(
            model_id=model_id, config_json=self.config_json
        )
        is_valid_url, _ = mdl_url_resolver.resolve_valid_hosted_model_url(model_id)
        return is_valid_url

    def is_available(self, model_id: str) -> bool:
        """
        Check if the model is available at the hosted URL.

        Parameters
        ----------
        model_id : str
            ID of the model to check.

        Returns
        -------
        bool
            True if the model is available, False otherwise.
        """
        if self.url is None:
            return self._is_available_unknown_url(model_id=model_id)
        else:
            return self._is_available_known_url()

    def _update_url(self, model_id: str):
        if self.url is None:
            from_hosted_file = os.path.join(
                self._model_path(model_id), IS_FETCHED_FROM_HOSTED_FILE
            )
            with open(from_hosted_file, "r") as f:
                data = json.load(f)
            self.url = data["url"]

    def write_apis(self, model_id: str):
        """
        Write APIs for the model.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        self.logger.debug("Writing APIs")
        di = HostedService(
            model_id=model_id, config_json=self.config_json, url=self.url
        )
        di.serve()
        di.close()

    def get_information(self, model_id: str):
        """
        Get information for the model.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        self.logger.debug(
            "Getting information for model identifier: {0}".format(model_id)
        )
        headers = {"accept": "*/*", "Content-Type": "application/json"}
        data = {}
        response = requests.post(
            self.url + "/info", headers=headers, data=json.dumps(data)
        )
        info = response.json()
        info_file = os.path.join(EOS, "dest", model_id, INFORMATION_FILE)
        with open(info_file, "w") as f:
            json.dump(info, f, indent=4)

    def get_metadata(self, model_id: str):
        """
        Get metadata for the model.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        self.logger.debug(
            "Getting api_schema for model identifier: {0}".format(model_id)
        )
        info_file = os.path.join(EOS, "dest", model_id, INFORMATION_FILE)
        with open(info_file, "r") as f:
            info = json.load(f)
        api_schema = info["api_schema"]
        api_schema_file = os.path.join(EOS, "dest", model_id, API_SCHEMA_FILE)
        with open(api_schema_file, "w") as f:
            json.dump(api_schema, f, indent=4)

    def write_status(self, model_id: str):
        """
        Write the status file for the model.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        status = {"done": True}
        status_file = os.path.join(EOS, "dest", model_id, STATUS_FILE)
        with open(status_file, "w") as f:
            json.dump(status, f, indent=4)

    def fetch(self, model_id: str):
        """
        Fetch the model from the hosted URL.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        self.logger.debug(
            "Fetching from hosted, model identifier: {0}".format(model_id)
        )
        mr = ModelRegisterer(model_id=model_id, config_json=self.config_json)
        self.logger.debug("Registering model")
        mr.register(is_from_hosted=True)
        self.logger.debug("Writing APIs")
        self._update_url(model_id)
        self.write_apis(model_id)
        self.get_information(model_id)
        self.get_metadata(model_id)
        self.write_status(model_id)
