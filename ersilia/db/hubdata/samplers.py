import os
import csv
import json
import requests
import random

from ersilia.utils.exceptions_utils.card_exceptions import InputBaseInformationError
from .interfaces import JsonModelsInterface
from ... import ErsiliaBase


from ...utils.paths import get_metadata_from_base_dir

_MODEL_STATUS_READY = "Ready"
_STATUS_FIELD = "Status"
_MODEL_ID_FIELD = "Identifier"
_INPUT_TYPE_FIELD = "Input"
_INPUT_SHAPE_FIELD = "Input Shape"

_ERSILIA_MAINTAINED_INPUTS_GITHUB_REPOSITORY = "ersilia-model-hub-maintained-inputs"


class ModelSampler(ErsiliaBase):
    """
    Samples random working models from the S3 model hub for downstream automations.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for initializing the sampler.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _get_models_from_s3_models_json(self):
        json_models_interface = JsonModelsInterface(config_json=self.config_json)
        models = json_models_interface.items_all()
        model_ids = [
            model[_MODEL_ID_FIELD]
            for model in models
            if model[_STATUS_FIELD] == _MODEL_STATUS_READY
        ]
        return model_ids

    def sample(self, n_samples, file_name=None):
        """
        Samples a specified number of models from S3 model hub.

        Parameters
        ----------
        n_samples : int
            Number of models to sample.
        file_name : str, optional
            Name of the file to save the sampled models.

        Returns
        -------
        list
            List of sampled model identifiers.
        """
        mdls = self._get_models_from_s3_models_json()
        sampled = random.sample(mdls, min(len(mdls), n_samples))
        return sampled
