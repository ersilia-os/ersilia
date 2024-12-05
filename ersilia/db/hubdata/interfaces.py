import importlib
import json

import boto3
import requests

from ... import ErsiliaBase
from ...default import AIRTABLE_MODEL_HUB_BASE_ID, AIRTABLE_MODEL_HUB_TABLE_NAME, ERSILIA_MODEL_HUB_S3_BUCKET, MODELS_JSON
from ...setup.requirements.pyairtable import PyAirtableRequirement






AIRTABLE_MAX_ROWS = 100000
AIRTABLE_PAGE_SIZE = 100


class AirtableInterface(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.base_id = AIRTABLE_MODEL_HUB_BASE_ID
        self.table_name = AIRTABLE_MODEL_HUB_TABLE_NAME
        self.max_rows = AIRTABLE_MAX_ROWS
        self.page_size = AIRTABLE_PAGE_SIZE
        self.write_api_key = None
        self.table = self._create_table(api_key=self._get_read_only_airtable_api_key())

    def _create_table(self, api_key):
        pyairtable_req = PyAirtableRequirement()
        if not pyairtable_req.is_installed():
            self.logger.debug("Installing PyAirTable from pip")
            pyairtable_req.install()
        pyairtable = importlib.import_module("pyairtable")
        return pyairtable.Table(api_key, self.base_id, self.table_name)

    @staticmethod
    def _get_read_only_airtable_api_key():
        url = "https://ersilia-model-hub.s3.eu-central-1.amazonaws.com/read_only_keys.json"
        r = requests.get(url)
        data = r.json()
        return data["AIRTABLE_READONLY_API_KEY"]

    def set_write_api_key(self, write_api_key):
        self.write_api_key = write_api_key
        self.table = self._create_table(api_key=write_api_key)

    def items(self):
        for records in self.table.iterate(
            page_size=self.page_size, max_records=self.max_rows
        ):
            for record in records:
                yield record

    def items_all(self):
        records = self.table.all()
        for record in records:
            yield record


class JsonModelsInterface(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        # self.cache_dir =
        self.json_file_name = MODELS_JSON
        self.url = f"https://{ERSILIA_MODEL_HUB_S3_BUCKET}.s3.eu-central-1.amazonaws.com/{MODELS_JSON}"

    def _read_json_file(self):
        response = requests.get(self.url)
        models_list = response.json()
        return models_list

    def items(self):
        models = self._read_json_file()
        for mdl in models:
            yield mdl

    def items_all(self):
        models = self._read_json_file()
        return models

    def identifier_exists(self, model_id):
        data = self._read_json_file()
        return any(item["Identifier"] == model_id for item in data)
