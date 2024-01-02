import requests
from pyairtable import Table
from ... import ErsiliaBase
from ...default import AIRTABLE_MODEL_HUB_BASE_ID, AIRTABLE_MODEL_HUB_TABLE_NAME

AIRTABLE_MAX_ROWS = 100000
AIRTABLE_PAGE_SIZE = 100


class AirtableInterface(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.api_key = self._get_read_only_airtable_api_key()
        self.base_id = AIRTABLE_MODEL_HUB_BASE_ID
        self.table_name = AIRTABLE_MODEL_HUB_TABLE_NAME
        self.max_rows = AIRTABLE_MAX_ROWS
        self.page_size = AIRTABLE_PAGE_SIZE
        self.write_api_key = None
        self.table = Table(self.api_key, self.base_id, self.table_name)

    @staticmethod
    def _get_read_only_airtable_api_key():
        url = "https://raw.githubusercontent.com/ersilia-os/ersilia/master/config/read_only_keys.json"
        r = requests.get(url)
        data = r.json()
        return data["AIRTABLE_READONLY_API_KEY"]

    def set_write_api_key(self, write_api_key):
        self.write_api_key = write_api_key
        self.table = Table(self.write_api_key, self.base_id, self.table_name)

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
