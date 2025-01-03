from ....serve.autoservice import AutoService
from . import BaseAction


class ModelChecker(BaseAction):
    """
    Checks the model by running an autoservice class.

    Parameters
    ----------
    model_id : str
        The model identifier.
    config_json : dict
        The configuration settings in JSON format.
    """

    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )

    def check(self):
        """
        Check that the autoservice works.
        """
        self.logger.debug("Checking that autoservice works")
        AutoService(self.model_id)
