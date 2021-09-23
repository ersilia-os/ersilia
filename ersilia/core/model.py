import os
import json
from ..core.modelbase import ModelBase
from ..serve.autoservice import AutoService
from ..serve.schema import ApiSchema
from ..io.input import ExampleGenerator
from ..default import MODEL_SIZE_FILE, CARD_FILE
from .. import logger


class ErsiliaModel(AutoService):
    def __init__(self, model, config_json=None, overwrite=False, verbose=False):
        self.logger = logger
        if verbose:
            self.logger.set_verbosity(1)
        else:
            self.logger.set_verbosity(0)
        model = ModelBase(model)
        self.overwrite = overwrite
        self.config_json = config_json
        self.model_id = model.model_id
        self.slug = model.slug
        self.api_schema = ApiSchema(model_id=self.model_id, config_json=self.config_json)
        AutoService.__init__(self, self.model_id, config_json=self.config_json)

    @property
    def input_type(self):
        with open(os.path.join(self._model_path(self.model_id), CARD_FILE), "r") as f:
            return [x.lower() for x in json.load(f)["Input"]]

    @property
    def output_type(self):
        with open(os.path.join(self._model_path(self.model_id), CARD_FILE), "r") as f:
            return [x.lower() for x in json.load(f)["Output"]]

    @property
    def schema(self):
        return self.api_schema.schema

    @property
    def size(self):
        with open(
            os.path.join(self._model_path(self.model_id), MODEL_SIZE_FILE), "r"
        ) as f:
            return json.load(f)

    def example(self, n_samples, file_name=None, simple=True):
        eg = ExampleGenerator(model_id=self.model_id, config_json=self.config_json)
        return eg.example(n_samples=n_samples, file_name=file_name, simple=simple)
