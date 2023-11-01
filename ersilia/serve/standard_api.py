import os
import csv
import json
import importlib
import requests

from .. import ErsiliaBase
from ..default import (
    EXAMPLE_STANDARD_INPUT_CSV_FILENAME,
    EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME,
)
from ..default import INFORMATION_FILE
from ..default import DEFAULT_API_NAME

MAX_INPUT_ROWS_STANDARD = 1000


class StandardCSVRunApi(ErsiliaBase):
    def __init__(self, model_id, url, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.logger.info(
            "You are running the app with a standard runner. Beware that this runner does not do as many checks on the input as the conventional runner: use it at your own risk."
        )
        self.model_id = model_id
        if url[-1] == "/":
            self.url = url[:-1]
        else:
            self.url = url
        self.logger.debug("Standard API processor started at {0}".format(self.url))
        self.api_name = DEFAULT_API_NAME
        self.path = os.path.abspath(self._model_path(self.model_id))
        self.standard_input_csv = os.path.join(
            self.path, EXAMPLE_STANDARD_INPUT_CSV_FILENAME
        )
        self.standard_output_csv = os.path.join(
            self.path, EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME
        )
        self.input_type = self.get_input_type()
        self.logger.debug("This is the input type: {0}".format(self.input_type))
        self.encoder = self.get_identifier_object_by_input_type()
        self.header = self.get_expected_output_header(self.standard_output_csv)
        self.logger.debug(
            "This is the expected header (max 10): {0}".format(self.header[:10])
        )

    def get_identifier_object_by_input_type(self):
        identifier_module_path = "ersilia.utils.identifiers.{0}".format(
            self.input_type[0].lower()
        )
        identifier_object = importlib.import_module(identifier_module_path).Identifier()
        return identifier_object

    def is_ready(self):
        if not self.is_input_type_standardizable():
            return False
        if not os.path.exists(self.standard_output_csv) or not os.path.exists(
            self.standard_input_csv
        ):
            return False
        else:
            return True

    def _is_input_file_too_long(self, input_data):
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            next(reader)
            c = 0
            for _ in reader:
                c += 1
                if c > MAX_INPUT_ROWS_STANDARD:
                    return True
        return False

    def is_input_standard_csv_file(self, input_data):
        if type(input_data) != str:
            return False
        if not input_data.endswith(".csv"):
            return False
        if not os.path.exists(input_data):
            return False
        if self._is_input_file_too_long(input_data):
            return False
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            header = next(reader)
        if len(header) == 1:
            h = header[0].lower()
            if not self.encoder.is_input_header(h):
                return False
            else:
                return True
        elif len(header) == 2:
            h0 = header[0].lower()
            h1 = header[1].lower()
            if not self.encoder.is_key_header(h0):
                return False
            if not self.encoder.is_input_header(h1):
                return False
            return True
        elif len(header) == 3:
            h0 = header[0].lower()
            h1 = header[1].lower()
            if not self.encoder.is_key_header(h0):
                return False
            if not self.encoder.is_input_header(h1):
                return False
            return True
        else:
            return False

    def get_input_type(self):
        with open(os.path.join(self.path, INFORMATION_FILE), "r") as f:
            info = json.load(f)
            return info["metadata"]["Input"]

    def is_input_type_standardizable(self):
        if len(self.input_type) != 1:
            return False
        return True

    def is_output_type_standardizable(self):
        with open(os.path.join(self.path, INFORMATION_FILE), "r") as f:
            api_schema = json.load(f)["api_schema"]
            if DEFAULT_API_NAME not in api_schema:
                return False
            meta = api_schema[DEFAULT_API_NAME]
            output_keys = meta["output"].keys()
            if len(output_keys) != 1:
                return False
        return True

    def is_output_csv_file(self, output_data):
        if type(output_data) != str:
            return False
        if not output_data.endswith(".csv"):
            return False
        return True

    def get_expected_output_header(self, output_data):
        with open(output_data, "r") as f:
            reader = csv.reader(f)
            header = next(reader)
        return header

    def serialize_to_json_three_columns(self, input_data):
        json_data = []
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for r in reader:
                json_data += [{"key": r[0], "input": r[1], "text": r[2]}]
        return json_data

    def serialize_to_json_two_columns(self, input_data):
        json_data = []
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for r in reader:
                json_data += [{"key": r[0], "input": r[1], "text": r[1]}]
        return json_data

    def serialize_to_json_one_column(self, input_data):
        json_data = []
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for r in reader:
                key = self.encoder.encode(r[0])
                json_data += [{"key": key, "input": r[0], "text": r[0]}]
        return json_data

    def serialize_to_json(self, input_data):
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            h = next(reader)
        if len(h) == 1:
            self.logger.debug("One column found in input")
            return self.serialize_to_json_one_column(input_data=input_data)
        elif len(h) == 2:
            self.logger.debug("Two columns found in input")
            return self.serialize_to_json_two_columns(input_data=input_data)
        elif len(h) == 3:
            self.logger.debug("Three columns found in input")
            return self.serialize_to_json_three_columns(input_data=input_data)
        else:
            self.logger.info(
                "More than two columns found in input! This is not standard."
            )
            return None

    def is_amenable(self, input_data, output_data):
        if not self.is_input_type_standardizable():
            return False
        if not self.is_output_type_standardizable():
            return False
        if not self.is_input_standard_csv_file(input_data):
            return False
        if not self.is_output_csv_file(output_data):
            return False
        self.logger.debug("It seems amenable for standard run")
        return True

    def serialize_to_csv(self, input_data, result, output_data):
        k = list(result[0].keys())[0]
        v = result[0][k]
        if type(v) is list:
            is_list = True
        else:
            is_list = False
        with open(output_data, "w") as f:
            writer = csv.writer(f)
            writer.writerow(self.header)
            for i_d, r_d in zip(input_data, result):
                v = r_d[k]
                if not is_list:
                    r = [i_d["key"], i_d["input"]] + [v]
                else:
                    r = [i_d["key"], i_d["input"]] + v
                writer.writerow(r)
        return output_data

    def post(self, input, output):
        input_data = self.serialize_to_json(input)
        url = "{0}/{1}".format(self.url, self.api_name)
        response = requests.post(url, json=input_data)
        if response.status_code == 200:
            result = response.json()
            output_data = self.serialize_to_csv(input_data, result, output)
            return output_data
        else:
            return None


class StandardQueryApi(object):
    def __init__(self, model_id, url):
        # TODO This class will be used to query directly the calculations lake.
        pass
