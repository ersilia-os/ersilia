import os
import shutil
import json
import datetime
import validators

from .... import ErsiliaBase
from .... import EOS
from .... import throw_ersilia_exception
from ....default import (
    IS_FETCHED_FROM_DOCKERHUB_FILE,
    IS_FETCHED_FROM_HOSTED_FILE,
    SERVICE_CLASS_FILE,
)
from ....db.hubdata.interfaces import AirtableInterface
from ....utils.exceptions_utils.fetch_exceptions import InvalidUrlError


class ModelRegisterer(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def register_from_dockerhub(self):
        data = {"docker_hub": True}
        self.logger.debug(
            "Registering model {0} in the file system".format(self.model_id)
        )
        path = os.path.join(EOS, "dest", self.model_id)
        self.logger.debug(path)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
        file_name = os.path.join(path, IS_FETCHED_FROM_DOCKERHUB_FILE)
        self.logger.debug(file_name)
        with open(file_name, "w") as f:
            json.dump(data, f)
        current_time = datetime.datetime.now()
        folder_name = current_time.strftime("%Y%m%d%H%M%S")
        path = os.path.join(EOS, "repository", self.model_id)
        if os.path.exists(path):
            shutil.rmtree(path)
        path = os.path.join(path, folder_name)
        os.makedirs(path)
        file_name = os.path.join(path, IS_FETCHED_FROM_DOCKERHUB_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)
        file_name = os.path.join(path, SERVICE_CLASS_FILE)
        self.logger.debug("Writing service class pulled_docker {0}".format(file_name))
        with open(file_name, "w") as f:
            f.write("pulled_docker")

    def register_not_from_dockerhub(self):
        data = {"docker_hub": False}
        path = self._model_path(self.model_id)
        file_name = os.path.join(path, IS_FETCHED_FROM_DOCKERHUB_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)
        path = self._get_bundle_location(model_id=self.model_id)
        file_name = os.path.join(path, IS_FETCHED_FROM_DOCKERHUB_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)

    def _resolve_url(self):
        self.logger.debug("Trying to find an available URL where the model is hosted")
        url_field = "Host URL"
        identifier_field = "Identifier"
        ai = AirtableInterface(config_json=self.config_json)
        for record in ai.items_all():
            fields = record["fields"]
            if fields[identifier_field] == self.model_id:
                if url_field not in fields:
                    self.logger.debug("No hosted URL found for this model")
                    return None
                url = fields[url_field]
                if validators.url(url):
                    self.logger.debug(
                        "This model has an associated URL: {0}".format(url)
                    )
                    return url
                else:
                    self.logger.debug(
                        "This doesn't seem to be a valid URL: {0}".format(url)
                    )
        self.logger.debug("Model was not found in AirTable")
        return None

    @throw_ersilia_exception
    def register_from_hosted(self, url=None):
        if url is None:
            url = self._resolve_url()
        else:
            if not validators.url(url):
                raise InvalidUrlError(url)
        data = {"hosted": True, "url": url}
        self.logger.debug(
            "Registering model {0} in the file system".format(self.model_id)
        )
        path = os.path.join(EOS, "dest", self.model_id)
        self.logger.debug(path)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
        file_name = os.path.join(path, IS_FETCHED_FROM_HOSTED_FILE)
        self.logger.debug(file_name)
        with open(file_name, "w") as f:
            json.dump(data, f)
        current_time = datetime.datetime.now()
        folder_name = current_time.strftime("%Y%m%d%H%M%S")
        path = os.path.join(EOS, "repository", self.model_id)
        if os.path.exists(path):
            shutil.rmtree(path)
        path = os.path.join(path, folder_name)
        os.makedirs(path)
        file_name = os.path.join(path, IS_FETCHED_FROM_HOSTED_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)
        file_name = os.path.join(path, SERVICE_CLASS_FILE)
        self.logger.debug("Writing service class hosted {0}".format(file_name))
        with open(file_name, "w") as f:
            f.write("hosted")

    def register_not_from_hosted(self):
        data = {"hosted": False}
        path = self._model_path(self.model_id)
        file_name = os.path.join(path, IS_FETCHED_FROM_HOSTED_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)
        path = self._get_bundle_location(model_id=self.model_id)
        file_name = os.path.join(path, IS_FETCHED_FROM_HOSTED_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)

    def register(self, is_from_dockerhub=False, is_from_hosted=False):
        if is_from_dockerhub and is_from_hosted:
            raise Exception
        if is_from_dockerhub and not is_from_hosted:
            self.register_from_dockerhub()
            self.register_not_from_dockerhub()
        if not is_from_dockerhub and is_from_hosted:
            self.register_from_hosted()
            self.register_not_from_dockerhub()
        else:
            self.register_not_from_dockerhub()
            self.register_not_from_hosted()
