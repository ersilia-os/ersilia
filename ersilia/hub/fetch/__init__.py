import validators
from ...db.hubdata.interfaces import AirtableInterface
from ...db.hubdata.json_models_interface import JsonModelsInterface
from ... import ErsiliaBase

MODEL_INSTALL_COMMANDS_FILE = "model_install_commands.sh"
DOCKERFILE = "Dockerfile"
ENVIRONMENT_YML = "environment.yml"
REQUIREMENTS_TXT = "requirements.txt"
STATUS_FILE = "status.json"
DONE_TAG = "done"

HOST_URL = "Host URL"
IDENTIFIER = "Identifier"


class ModelURLResolver(ErsiliaBase):
    def __init__(self, model_id, config_json=None, credentials_json=None):
        super().__init__(config_json, credentials_json)
        self.model_id = model_id
        self.ai = AirtableInterface(config_json=self.config_json)
        self.ji = JsonModelsInterface(config_json=self.config_json)

    def _find_url_using_s3_models_json(self, model_id):
        self.logger.debug(
            "Trying to find an available URL where the model is hosted using S3 Models JSON"
        )
        for mdl in self.ji.items():
            if mdl[IDENTIFIER] == model_id:
                if HOST_URL in mdl:
                    url = mdl[HOST_URL]
                    return url
                else:
                    self.logger.debug(
                        "No hosted URL found for this model in S3 Models JSON"
                    )
                    return None
        self.logger.debug("Model was not found in S3 Models JSON")

    def _find_url_using_airtable(self, model_id):
        url_field = HOST_URL
        identifier_field = IDENTIFIER
        for record in self.ai.items_all():
            fields = record["fields"]
            if fields[identifier_field] == model_id:
                if url_field in fields:
                    return fields[url_field]
                else:
                    self.logger.debug("No hosted URL found for this model in AirTable")
                    return
        self.logger.debug("Model was not found in AirTable")

    def resolve_valid_hosted_model_url(self, model_id):
        """Resolves the URL of a hosted model

        Parameters
        ----------
        model_id : str
            Alpha numeric identifier literal of the model

        Returns
        -------
        Tuple[bool, str]
            Tuple containing a boolean indicating if the URL is valid and the URL itself
        """
        url = self._find_url_using_s3_models_json(
            model_id
        ) or self._find_url_using_airtable(model_id)
        if url:
            if validators.url(url):
                self.logger.debug("This model has an associated URL: {0}".format(url))
                return (True, url)
            else:
                self.logger.debug(
                    "This doesn't seem to be a valid URL: {0}".format(url)
                )
        return (False, url)
