import os
import csv
import json
import importlib
import requests
import asyncio
import nest_asyncio
from ..store.api import InferenceStoreApi
from ..store.utils import OutputSource
from .. import ErsiliaBase
from ..default import (
    EXAMPLE_STANDARD_INPUT_CSV_FILENAME,
    EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME,
)
from ..default import INFORMATION_FILE, API_SCHEMA_FILE
from ..default import DEFAULT_API_NAME

MAX_INPUT_ROWS_STANDARD = 1000

nest_asyncio.apply()

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
        self.validate_smiles = self.get_identifier_object_by_input_type().validate_smiles
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
        try:
            with open(os.path.join(self.path, INFORMATION_FILE), "r") as f:
                info = json.load(f)
                if "metadata" in info and "Input" in info["metadata"]:
                    return info["metadata"]["Input"]
                elif "card" in info and "Input" in info["card"]:
                    return info["card"]["Input"]
                else:
                    raise KeyError("Neither 'metadata' nor 'card' contains 'Input' key.")
        
        except FileNotFoundError:
            self.logger.debug(f"Error: File '{INFORMATION_FILE}' not found in the path '{self.path}'")
        except json.JSONDecodeError:
            self.logger.debug(f"Error: Failed to parse JSON in file '{INFORMATION_FILE}'")
        except KeyError as e:
             self.logger.debug(f"Error: {e}")
        except Exception as e:
             self.logger.debug(f"An unexpected error occurred: {e}")


    def is_input_type_standardizable(self):
        if len(self.input_type) != 1:
            return False
        return True

    def is_output_type_standardizable(self):
        api_schema_file_path = os.path.join(self.path, API_SCHEMA_FILE)
        try:
            with open(api_schema_file_path, "r") as f:
                api_schema = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return False

        meta = api_schema.get(DEFAULT_API_NAME)
        if not meta or len(meta.get("output", {})) != 1:
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
     

    def parse_smiles_list(self, input_data):
        if not input_data or all(not s.strip() for s in input_data):
            raise ValueError("The list of SMILES strings is empty or contains only empty strings.")
        return [{'key': self.encoder.encode(smiles), 'input': smiles, 'text': smiles} for smiles in input_data if self.validate_smiles(smiles)]

    def parse_smiles_string(self, input):
        if not self.validate_smiles(input):
            raise ValueError("The SMILES string is invalid.")
        key = self.encoder.encode(input)
        return [{'key': key, 'input': input, 'text': input}]
    
    def serialize_to_json_three_columns(self, input_data):
        json_data = []
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            next(reader)  
            for row in reader:
                if self.validate_smiles(row[1]):  
                    json_data += [{"key": row[0], "input": row[1], "text": row[2]}]
        return json_data

    def serialize_to_json_two_columns(self, input_data):
        json_data = []
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                if self.validate_smiles(row[1]):  
                    json_data += [{"key": row[0], "input": row[1], "text": row[1]}]
        return json_data

    def serialize_to_json_one_column(self, input_data):
        json_data = []
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                if self.validate_smiles(row[0]):
                    key = self.encoder.encode(row[0])
                    json_data += [{"key": key, "input": row[0], "text": row[0]}]
        return json_data

    async def async_serialize_to_json_one_column(self, input_data):
        smiles_list = self.get_list_from_csv(input_data)
        smiles_list = [smiles for smiles in smiles_list if self.validate_smiles(smiles)]
        json_data = await self.encoder.encode_batch(smiles_list)
        return json_data

    def get_list_from_csv(self, input_data):
        smiles_list = []
        with open(input_data, mode='r') as file:
            reader = csv.DictReader(file)
            header = reader.fieldnames
            key = header[0] if len(header) == 1 else header[1]
            for row in reader:
                smiles = row.get(key)
                if smiles and smiles not in smiles_list and self.validate_smiles(smiles):
                    smiles_list.append(smiles)
        return smiles_list

    def serialize_to_json(self, input_data):
        if isinstance(input_data, str) and os.path.isfile(input_data):
            with open(input_data, "r") as f:
                reader = csv.reader(f)
                h = next(reader)
            if len(h) == 1:
                self.logger.debug("One column found in input")
                return asyncio.run(self.async_serialize_to_json_one_column(input_data))
            elif len(h) == 2:
                self.logger.debug("Two columns found in input")
                return self.serialize_to_json_two_columns(input_data=input_data)
            elif len(h) == 3:
                self.logger.debug("Three columns found in input")
                return self.serialize_to_json_three_columns(input_data=input_data)
            else:
                self.logger.info("More than three columns found in input! This is not standard.")
                return None
        elif isinstance(input_data, str):
            return self.parse_smiles_string(input_data)
        elif isinstance(input_data, list):
            return self.parse_smiles_list(input_data)
        else:
            raise ValueError("Input must be either a file path (string), a SMILES string, or a list of SMILES strings.")
    
    def is_amenable(self, output_data):
        if not self.is_input_type_standardizable():
            return False
        if not self.is_output_type_standardizable():
            return False
        if not self.is_output_csv_file(output_data):
            return False
        self.logger.debug("It seems amenable for standard run")
        return True
    
    def serialize_to_csv(self, input_data, result, output_data):
        with open(output_data, "w") as f:
            writer = csv.writer(f)
            writer.writerow(self.header)
            for i_d, r_d in zip(input_data, result):
                r = [i_d["key"], i_d["input"]]
                for k in self.header[2:]:
                    v = r_d[k]
                    if isinstance(v, list):
                        r+=v
                    else:
                        r+=[v]
                writer.writerow(r)
        return output_data

    def post(self, input, output, output_source=OutputSource.LOCAL_ONLY):
        input_data = self.serialize_to_json(input)
        if OutputSource.is_cloud(output_source):
            store = InferenceStoreApi(model_id=self.model_id)
            return store.get_precalculations(input_data)
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
