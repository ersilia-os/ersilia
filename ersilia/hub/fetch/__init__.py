import validators

from ... import ErsiliaBase
from ...db.hubdata.interfaces import JsonModelsInterface

MODEL_INSTALL_COMMANDS_FILE = "model_install_commands.sh"
DOCKERFILE = "Dockerfile"
ENVIRONMENT_YML = "environment.yml"
REQUIREMENTS_TXT = "requirements.txt"
STATUS_FILE = "status.json"
DONE_TAG = "done"

HOST_URL = "Host URL"
IDENTIFIER = "Identifier"


class ModelURLResolver(ErsiliaBase):
    """
    Class to resolve the URL of a model.

    This class provides methods to resolve the URL of a model based on its ID.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    config_json : dict, optional
        Configuration settings in JSON format.
    credentials_json : dict, optional
        Credentials settings in JSON format.
    """

    def __init__(self, model_id, config_json=None, credentials_json=None):
        super().__init__(config_json, credentials_json)
        self.model_id = model_id
        self.ji = JsonModelsInterface(config_json=self.config_json)
        self.models_cache = None

    def _cache_models(self):
        if self.models_cache is None:
            models = self.ji._read_json_file()
            self.models_cache = {mdl["Identifier"]: mdl for mdl in models}

    def _find_url_using_s3_models_json(self, model_id):
        self.logger.debug(
            "Trying to find an available URL where the model is hosted using S3 Models JSON"
        )
        self._cache_models()
        model = self.models_cache.get(model_id)

        if model:
            if "Host URL" in model:
                return model["Host URL"]
            else:
                self.logger.debug(
                    "No hosted URL found for this model in S3 Models JSON"
                )
                return None
        self.logger.debug("Model was not found in S3 Models JSON")

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
        url = self._find_url_using_s3_models_json(model_id)
        if url:
            if validators.url(url):
                self.logger.debug("This model has an associated URL: {0}".format(url))
                return (True, url)
            else:
                self.logger.debug(
                    "This doesn't seem to be a valid URL: {0}".format(url)
                )
        return (False, url)
