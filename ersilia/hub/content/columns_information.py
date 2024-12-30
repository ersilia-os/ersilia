import os
import csv
from urllib.request import urlopen
from ... import ErsiliaBase


ROOT = os.path.dirname(os.path.abspath(__file__))

EXPECTED_HEADER = [
    "name",
    "type",
    "direction",
    "description"
]

MIN_DESCRIPTION_LENGTH = 60


class ColumnsInformation(ErsiliaBase):
    def __init__(self, model_id, api_name, config_json=None):
        self.model_id = model_id
        self.api_name = api_name
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.relative_path = "model/framework/columns/{0}_columns.csv".format(self.api_name)
        with open(os.path.join(ROOT, "columns", "data_types.txt"), "r") as f:
            self.DATA_TYPES = []
            for l in f:
                self.DATA_TYPES += [l.strip()]
        with open(os.path.join(ROOT, "columns", "desired_directions.txt"), "r") as f:
            self.DESIRED_DIRECTIONS = []
            for l in f:
                self.DESIRED_DIRECTIONS += [l.strip()]

    def _get_columns_information_from_github(self):
        org = "ersilia-os"
        branch = "main"
        url = "https://raw.githubusercontent.org/{0}/{1}/{2}/{3}".format(org, self.model_id, branch, self.relative_path)
        # TODO try / except
        with urlopen(url) as response:
            pass
            # TODO read text

    def _get_columns_information_from_local(self):
        file_name = os.path.join(self._model_path(self.model_id), self.relative_path)
        if os.path.exists(file_name):
            with open(file_name, "r") as f:
                names = []
                types = []
                directions = []
                descriptions = []
                reader = csv.reader(f)
                header = next(reader)
                if header != EXPECTED_HEADER:
                    raise ValueError("Header {0} is not {1}".format(header, EXPECTED_HEADER))
                for r in reader:
                    names += [r[0]]
                    types += [r[1]]
                    if r[2] == "":
                        directions += [None]
                    else:
                        directions += [r[2]]
                    descriptions += [r[3]]
            return {"name": names, "description": descriptions, "direction": directions}
        else:
            self.logger.debug("Explicit columns data for {0} API does not exist in file {1}".format(self.api_name, file_name))
            return None

    def _validate_columns_data(self, data):
        for d in data["name"]:
            if d[0].lower() != d[0]:
                raise ValueError("Column names must be lowercase")
            if not d.replace("_", "").isalnum():
                raise ValueError("Column names must be alphanumeric or contain underscores")
        for d in data["description"]:
            if len(d) < MIN_DESCRIPTION_LENGTH:
                raise ValueError("Description is too short. A minimum of {0} characters is expected".format(MIN_DESCRIPTION_LENGTH))
        for d in data["direction"]:
            if d not in self.DESIRED_DIRECTIONS:
                raise ValueError("Direction {0} is not an accepted direction: {1}".format(d, self.DESIRED_DIRECTIONS))

    def load(self):
        data = self._get_columns_information_from_local()
        if data is None:
            data = self._get_columns_information_from_github()
        if data is None:
            return None
        self._validate_columns_data(data)
        return data
        
