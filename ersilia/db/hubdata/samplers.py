import os
import csv
import json
import requests
import random

from ersilia.utils.exceptions_utils.card_exceptions import InputBaseInformationError
from .interfaces import AirtableInterface
from .json_models_interface import JsonModelsInterface
from ... import ErsiliaBase


from ...utils.paths import get_metadata_from_base_dir

_MODEL_STATUS_READY = "Ready"
_STATUS_FIELD = "Status"
_MODEL_ID_FIELD = "Identifier"
_INPUT_TYPE_FIELD = "Input"
_INPUT_SHAPE_FIELD = "Input Shape"

_ERSILIA_MAINTAINED_INPUTS_GITHUB_REPOSITORY = "ersilia-model-hub-maintained-inputs"


class ModelSampler(ErsiliaBase):
    """Get a random working model from the model hub to use in downstream automations."""

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _get_models_from_airtable(self):
        airtable_interface = AirtableInterface(config_json=self.config_json)
        model_ids = []
        for record in airtable_interface.items():
            status = record["fields"][_STATUS_FIELD]
            if status == _MODEL_STATUS_READY:
                model_ids += [record["fields"][_MODEL_ID_FIELD]]
        return model_ids

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
        entities = self._get_models_from_s3_models_json()
        if len(entities) == 0:
            entities = self._get_models_from_airtable()
        entities = random.sample(entities, min(len(entities), n_samples))
        if file_name is None:
            for e in entities:
                print(e)  # TODO Change print to click?
        else:
            with open(file_name, "w") as f:
                writer = csv.writer(f)
                for r in entities:
                    writer.writerow([r])


class InputSampler(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        result = self._get_input_type_and_shape()
        assert result is not None
        self.input_type = result[0]
        self.input_shape = result[1]

    def _get_input_type_and_shape_from_airtable(self):
        airtable_interface = AirtableInterface(config_json=self.config_json)
        for records in airtable_interface.items():
            fields = records["fields"]
            model_id = fields[_MODEL_ID_FIELD]
            input_type = fields[_INPUT_TYPE_FIELD]
            input_shape = fields[_INPUT_SHAPE_FIELD]
            status = fields[_STATUS_FIELD]
            if status == _MODEL_STATUS_READY:
                if model_id == self.model_id:
                    return input_type, input_shape
        return None

    def _get_input_type_and_shape_from_metadata(self):
        dest_path = self._model_path(self.model_id)
        try:
            data = get_metadata_from_base_dir(dest_path)
        except FileNotFoundError:
            return None
        input_type = data[_INPUT_TYPE_FIELD]
        input_shape = data[_INPUT_SHAPE_FIELD]
        return input_type, input_shape

    def _get_input_type_and_shape_from_s3_models_json(self):
        json_models_interface = JsonModelsInterface(config_json=self.config_json)
        for mdl in json_models_interface.items():
            model_id = mdl[_MODEL_ID_FIELD]
            input_type = mdl[_INPUT_TYPE_FIELD]
            input_shape = mdl[_INPUT_SHAPE_FIELD]
            status = mdl[_STATUS_FIELD]
            if status == _MODEL_STATUS_READY:
                if model_id == self.model_id:
                    return input_type, input_shape
        return None

    def _get_input_type_and_shape(self):
        res = self._get_input_type_and_shape_from_metadata()
        if res is not None:
            return res
        res = self._get_input_type_and_shape_from_s3_models_json()
        if res is not None:
            return res
        res = self._get_input_type_and_shape_from_airtable()
        if res is not None:
            return res
        return None

    def _create_url_to_get_sample_content(self):
        shapes = ["single", "pair", "list", "pair-of-lists"]
        # Ensure that there is only one input type specified in the metadata.
        assert len(self.input_type) == 1
        input_type = self.input_type[0].lower()
        input_shape = self.input_shape.lower().replace(" ", "-")
        if input_shape not in shapes:
            raise InputBaseInformationError()
        return "https://raw.githubusercontent.com/ersilia-os/{0}/main/{1}/{2}/inp-000.csv".format(
            _ERSILIA_MAINTAINED_INPUTS_GITHUB_REPOSITORY,
            input_type,
            input_shape,
        )

    def _get_inputs_from_maintained_file(self):
        url = self._create_url_to_get_sample_content()
        with requests.Session() as s:
            download = s.get(url)
            decoded_content = download.content.decode("utf-8")
            cr = csv.reader(decoded_content.splitlines(), delimiter=",")
            data = list(cr)
            headers = data[0]
            content = data[1:]
        return content, headers

    def sample(self, n_samples, file_name):
        entities, headers = self._get_inputs_from_maintained_file()
        entities = random.sample(entities, min(n_samples, len(entities)))
        if file_name is None:
            for e in entities:
                print(e[1])  # TODO Change print to click?
        else:
            with open(file_name, "w") as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                for r in entities:
                    writer.writerow(r)
