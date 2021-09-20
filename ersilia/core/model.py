from ..core.modelbase import ModelBase
from ..serve.autoservice import AutoService
from .. import logger


class ErsiliaModel(AutoService):
    def __init__(self, model, config_json=None, overwrite=False, verbose=False):
        if verbose:
            logger.set_verbosity(1)
        else:
            logger.set_verbosity(0)
        model = ModelBase(model)
        self.overwrite = overwrite
        self.config_json = config_json
        self.model_id = model.model_id
        self.slug = model.slug
        AutoService.__init__(self, self.model_id, config_json=config_json)
