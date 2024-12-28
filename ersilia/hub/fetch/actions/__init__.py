from ....core.base import ErsiliaBase


class BaseAction(ErsiliaBase):
    """
    Base class for actions.

    This class provides common methods for actions.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    config_json : dict
        Configuration settings in JSON format.
    credentials_json : dict
        Credentials settings in JSON format.
    """

    def __init__(self, model_id, config_json, credentials_json):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.model_id = model_id
