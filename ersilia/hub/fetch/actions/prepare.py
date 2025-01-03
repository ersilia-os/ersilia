from .... import throw_ersilia_exception
from ....utils.exceptions_utils.delete_exceptions import ModelDeleteError
from ...bundle.status import ModelStatus
from ...delete.delete import ModelFullDeleter
from . import BaseAction


class ModelPreparer(BaseAction):
    """
    Prepares a model for use by deleting existing data if necessary.

    Parameters
    ----------
    model_id : str
        Identifier of the model to be prepared.
    overwrite : bool
        Whether to overwrite existing data.
    config_json : dict
        Configuration settings for the preparer.

    Methods
    -------
    prepare()
        Prepares the model by deleting existing data if necessary.
    """

    def __init__(self, model_id: str, overwrite: bool, config_json: dict):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.overwrite = overwrite
        self.status = ModelStatus(config_json=self.config_json)
        self.deleter = ModelFullDeleter(
            config_json=self.config_json, overwrite=self.overwrite
        )

    @throw_ersilia_exception()
    def prepare(self):
        """
        Prepares the model by deleting existing data if necessary.
        """
        try:
            self.deleter.delete(self.model_id)
        except:
            raise ModelDeleteError(model=self.model_id)
