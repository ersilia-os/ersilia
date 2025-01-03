import asyncio
import csv
import importlib
import json
import os

import nest_asyncio
import requests

from .. import ErsiliaBase
from ..default import (
    API_SCHEMA_FILE,
    DEFAULT_API_NAME,
    EXAMPLE_STANDARD_INPUT_CSV_FILENAME,
    EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME,
    INFORMATION_FILE,
    PREDEFINED_EXAMPLE_FILES,
)
from ..store.api import InferenceStoreApi
from ..store.utils import OutputSource

MAX_INPUT_ROWS_STANDARD = 1000

nest_asyncio.apply()


class StandardCSVRunApi(ErsiliaBase):
    """
    Class for running standard CSV API. An API that performs few checks on the input data.

    Parameters
    ----------
    model_id : str
        The ID of the model to be used.
    url : str
        The URL of the API.
    config_json : dict, optional
        Configuration settings in JSON format.

    Examples
    --------
    .. code-block:: python
        model_id = "eosxxxx"
        url = "http://0.0.0.0:15221/run"
        api = StandardCSVRunApi(model_id, url)
        input_data = "path/to/input.csv"
        output_data = "path/to/output.csv"
        result = api.post(input_data, output_data)
    """

    def __init__(self, model_id, url, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.logger.info(
            "You are running the app with a standard runner. Beware that this runner does not do as many checks on the input as the conventional runner: use it at your own risk."
        )
        self.model_id = model_id
        # TODO WHY? Why can't we init with a cleaner url?
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
        metadata = self._read_information_file()
        self.input_type = self._read_field_from_metadata(metadata, "Input")
        self.input_shape = self._read_field_from_metadata(metadata, "Input Shape")
        self.logger.debug("This is the input type: {0}".format(self.input_type))
        self.encoder = self.get_identifier_object_by_input_type()
        # TODO This whole validate_smiles thing can go away since we already handle this in the encoder
        self.validate_smiles = (
            self.get_identifier_object_by_input_type().validate_smiles
        )  # TODO this can just be self.encoder.validate_smiles
        self.header = self.get_expected_output_header()
        self.logger.debug(
            "This is the expected header (max 10): {0}".format(self.header[:10])
        )

    def _read_information_file(self):
        try:
            with open(os.path.join(self.path, INFORMATION_FILE), "r") as f:
                info = json.load(f)
                return info
        except FileNotFoundError:
            self.logger.debug(
                f"Error: File '{INFORMATION_FILE}' not found in the path '{self.path}'"
            )
        except json.JSONDecodeError:
            self.logger.debug(
                f"Error: Failed to parse JSON in file '{INFORMATION_FILE}'"
            )
        except Exception as e:
            self.logger.debug(f"An unexpected error occurred: {e}")

    def _read_field_from_metadata(self, meta, field):
        if not meta:
            self.logger.error("No metadata given")
            return None
        if "metadata" in meta and field in meta["metadata"]:
            return meta["metadata"][field]
        elif "card" in meta and field in meta["card"]:
            return meta["card"][field]
        else:
            self.logger.error(f"Neither 'metadata' nor 'card' contains '{field}' key.")

    def get_identifier_object_by_input_type(self):
        """
        Get the identifier object by input type.

        Returns
        -------
        object
            The identifier object.
        """
        identifier_module_path = "ersilia.utils.identifiers.{0}".format(
            self.input_type[0].lower()
        )
        identifier_object = importlib.import_module(identifier_module_path).Identifier()
        return identifier_object

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

    def is_input_type_standardizable(self):
        """
        Check if the input type is standardizable.

        Returns
        -------
        bool
            True if the input type is standardizable, False otherwise.
        """
        if self.input_type and self.input_shape:
            if self.input_type[0] == "Compound" and self.input_shape == "Single":
                return True
        return False

    def is_output_type_standardizable(self):
        """
        Check if the output type is standardizable.

        Returns
        -------
        bool
            True if the output type is standardizable, False otherwise.
        """
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
        """
        Check if the output data is a CSV file.

        Parameters
        ----------
        output_data : str
            The output data to check.

        Returns
        -------
        bool
            True if the output data is a CSV file, False otherwise.
        """
        if type(output_data) != str:
            return False
        if not output_data.endswith(".csv"):
            return False
        return True

    def get_expected_output_header(self):
        """
        Calculate the expected output header from the predefined example files or the standard output file.

        Returns
        -------
        List | None
           Returns the header which is a list of column names, or None if the header could not be determined.
        """
        file = None
        for pf in PREDEFINED_EXAMPLE_FILES:
            if os.path.exists(os.path.join(self.path, pf)):
                file = os.path.join(self.path, pf)
                self.logger.debug(
                    f"Determining header from predefined example file: {pf}"
                )
                break
        if not file and os.path.exists(self.standard_output_csv):
            file = self.standard_output_csv
            self.logger.debug(
                f"Determining header from standard output file: {self.standard_output_csv}"
            )

        try:
            with open(file, "r") as f:
                reader = csv.reader(f)
                header = next(reader)
                if (
                    header[0:2]
                    != [
                        "key",
                        "input",
                    ]
                ):  # Slicing doesn't raise an error even if the list does not have 2 elements
                    header = ["key", "input"] + header
            return header
        except (FileNotFoundError, StopIteration):
            self.logger.error(f"Could not determine header from file {file}")

    def parse_smiles_list(self, input_data):
        """
        Parse a list of SMILES strings.

        Parameters
        ----------
        input_data : list
            List of SMILES strings.

        Returns
        -------
        list
            List of dictionaries containing encoded SMILES strings.
        """
        if not input_data or all(not s.strip() for s in input_data):
            raise ValueError(
                "The list of SMILES strings is empty or contains only empty strings."
            )
        return [
            {"key": self.encoder.encode(smiles), "input": smiles, "text": smiles}
            for smiles in input_data
            if self.validate_smiles(smiles)
        ]

    def parse_smiles_string(self, input):
        """
        Parse a single SMILES string.

        Parameters
        ----------
        input : str
            A SMILES string.

        Returns
        -------
        list
            List containing a dictionary with the encoded SMILES string.
        """
        if not self.validate_smiles(input):
            raise ValueError("The SMILES string is invalid.")
        key = self.encoder.encode(input)
        return [{"key": key, "input": input, "text": input}]

    def serialize_to_json_three_columns(self, input_data):
        """
        Serialize data to JSON with three columns.

        Parameters
        ----------
        input_data : str
            The input data file path.

        Returns
        -------
        list
            The serialized JSON data.
        """
        json_data = []
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                if self.validate_smiles(row[1]):
                    json_data += [{"key": row[0], "input": row[1], "text": row[2]}]
        return json_data

    def serialize_to_json_two_columns(self, input_data):
        """
        Serialize data to JSON with two columns.

        Parameters
        ----------
        input_data : str
            The input data file path.

        Returns
        -------
        list
            The serialized JSON data.
        """
        json_data = []
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                if self.validate_smiles(row[1]):
                    json_data += [{"key": row[0], "input": row[1], "text": row[1]}]
        return json_data

    def serialize_to_json_one_column(self, input_data):
        """
        Serialize data to JSON with one column.

        Parameters
        ----------
        input_data : str
            The input data file path.

        Returns
        -------
        list
            The serialized JSON data.
        """
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
        """
        Asynchronously serialize data to JSON with one column.

        Parameters
        ----------
        input_data : str
            The input data file path.

        Returns
        -------
        list
            The serialized JSON data.
        """
        smiles_list = self.get_list_from_csv(input_data)
        smiles_list = [smiles for smiles in smiles_list if self.validate_smiles(smiles)]
        json_data = await self.encoder.encode_batch(smiles_list)
        return json_data

    def get_list_from_csv(self, input_data):
        """
        Get a list from a CSV file.

        Parameters
        ----------
        input_data : str
            The input data file path.

        Returns
        -------
        list
            The list of data from the CSV file.
        """
        smiles_list = []
        with open(input_data, mode="r") as file:
            reader = csv.DictReader(file)
            header = reader.fieldnames
            key = header[0] if len(header) == 1 else header[1]
            for row in reader:
                smiles = row.get(key)
                if (
                    smiles
                    and smiles not in smiles_list
                    and self.validate_smiles(smiles)
                ):
                    smiles_list.append(smiles)
        return smiles_list

    def serialize_to_json(self, input_data):
        """
        Serialize input data to JSON format.

        Parameters
        ----------
        input_data : str | list
            Input data which can be a file path, a SMILES string, or a list of SMILES strings.

        Returns
        -------
        list
            List of dictionaries containing serialized input data.
        """
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
                self.logger.info(
                    "More than three columns found in input! This is not standard."
                )
                return None
        elif isinstance(input_data, str):
            return self.parse_smiles_string(input_data)
        elif isinstance(input_data, list):
            return self.parse_smiles_list(input_data)
        else:
            raise ValueError(
                "Input must be either a file path (string), a SMILES string, or a list of SMILES strings."
            )

    def is_amenable(self, output_data):
        """
        Check if the output data is amenable for a standard run.

        Parameters
        ----------
        output_data : str
            The output data to check.

        Returns
        -------
        bool
            True if the output data is amenable for a standard run, False otherwise.
        """
        if not self.header:
            self.logger.debug("Not amenable for standard run: header not found")
            return False
        if not self.is_input_type_standardizable():
            self.logger.debug(
                "Not amenable for standard run: input type not standardizable"
            )
            return False
        if not self.is_output_type_standardizable():
            self.logger.debug(
                "Not amenable for standard run: output type not standardizable"
            )
            return False
        if not self.is_output_csv_file(output_data):
            self.logger.debug("Not amenable for standard run: output data not CSV file")
            return False
        self.logger.debug("It seems amenable for standard run")
        return True

    def serialize_to_csv(self, input_data, result, output_data):
        """
        Serialize input data and result to a CSV file.

        Parameters
        ----------
        input_data : list
            List of dictionaries containing input data.
        result : list
            List of dictionaries containing result data.
        output_data : str
            Path to the output CSV file.

        Returns
        -------
        str
            Path to the output CSV file.
        """
        if isinstance(result, dict) and not list(result.keys()) == self.header:
            result = result[list(result.keys())[0]]
            if (
                isinstance(result[0], dict)
                and not list(result[0].keys()) == self.header
            ):
                for idx, item in enumerate(result):
                    result[idx] = item[list(item.keys())[0]]

        assert len(input_data) == len(result)

        with open(output_data, "w") as f:
            writer = csv.writer(f)
            writer.writerow(self.header)
            for i_d, r_d in zip(input_data, result):
                r = [i_d["key"], i_d["input"]]
                for k in self.header[2:]:
                    if isinstance(r_d, dict):
                        v = r_d[k]
                    else:
                        v = r_d
                    if isinstance(v, list):
                        r += v
                    else:
                        r += [v]
                writer.writerow(r)
        return output_data

    def post(self, input, output, output_source=OutputSource.LOCAL_ONLY):
        """
        Post input data to the API and get the output data. And store the output data if needed.

        Parameters
        ----------
        input : str | list
            Input data which can be a file path, a SMILES string, or a list of SMILES strings.
        output : str
            Path to the output CSV file.
        output_source : OutputSource, optional
            Source of the output data. Default is OutputSource.LOCAL_ONLY.

        Returns
        -------
        str | None
            Path to the output CSV file if successful, None otherwise.
        """
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
