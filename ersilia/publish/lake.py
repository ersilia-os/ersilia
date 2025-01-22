from .. import ErsiliaBase
from ..lake.manager import IsauraManager


class LakeStorer(ErsiliaBase):
    """
    Class to handle storing data in the lake.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    config_json : dict
        Configuration in JSON format.
    credentials_json : dict
        Credentials in JSON format.
    """

    def __init__(self, model_id, config_json, credentials_json):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.model_id = model_id
        self.isaura_manager = IsauraManager(
            model_id=model_id,
            config_json=config_json,
            credentials_json=credentials_json,
        )

    def store(self):
        """
        Store data in the lake.
        """
        self.logger.debug("Appeding local to public")
        self.isaura_manager.append_local_to_public()
        self.logger.debug("Pushing")
        self.isaura_manager.push()
