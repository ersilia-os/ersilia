import os
import requests
import json
import validators

from ..register.register import ModelRegisterer
from ....serve.services import HostedService
from ....db.hubdata.interfaces import AirtableInterface

from .... import ErsiliaBase
from .... import EOS
from ....default import API_SCHEMA_FILE, INFORMATION_FILE, IS_FETCHED_FROM_HOSTED_FILE
from .. import STATUS_FILE


class ModelHostedFetcher(ErsiliaBase):
    def __init__(self, url, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.logger.debug("Initialized with URL: {0}".format(url))
        self.url = url

    def _is_available_known_url(self):
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

    def _is_available_unknown_url(self, model_id):
        self.logger.debug("Trying to find an available URL where the model is hosted")
        url_field = "Host URL"
        identifier_field = "Identifier"
        ai = AirtableInterface(config_json=self.config_json)
        for record in ai.items_all():
            fields = record["fields"]
            if fields[identifier_field] == model_id:
                if url_field not in fields:
                    self.logger.debug("No hosted URL found for this model")
                    return False
                url = fields[url_field]
                if validators.url(url):
                    self.logger.debug(
                        "This model has an associated URL: {0}".format(url)
                    )
                    return True
                else:
                    self.logger.debug(
                        "This doesn't seem to be a valid URL: {0}".format(url)
                    )
        self.logger.debug("Model was not found in AirTable")
        return False

    def is_available(self, model_id):
        if self.url is None:
            return self._is_available_unknown_url(model_id=model_id)
        else:
            return self._is_available_known_url()

    def _update_url(self, model_id):
        if self.url is None:
            from_hosted_file = os.path.join(
                self._model_path(model_id), IS_FETCHED_FROM_HOSTED_FILE
            )
            with open(from_hosted_file, "r") as f:
                data = json.load(f)
            self.url = data["url"]

    def write_apis(self, model_id):
        self.logger.debug("Writing APIs")
        di = HostedService(
            model_id=model_id, config_json=self.config_json, url=self.url
        )
        di.serve()
        di.close()

    def get_information(self, model_id):
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

    def get_metadata(self, model_id):
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

    def write_status(self, model_id):
        status = {"done": True}
        status_file = os.path.join(EOS, "dest", model_id, STATUS_FILE)
        with open(status_file, "w") as f:
            json.dump(status, f, indent=4)

    def fetch(self, model_id):
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
