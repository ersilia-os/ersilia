import asyncio
import csv
import importlib
import json
import os
import time

import nest_asyncio
import requests

from .. import ErsiliaBase
from ..default import (
    DEFAULT_API_NAME,
    EXAMPLE_STANDARD_INPUT_CSV_FILENAME,
    EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME,
    INFORMATION_FILE,
)
from ..hub.content.columns_information import ColumnsInformation
from ..io.output import GenericOutputAdapter
from ..store.api import InferenceStoreApi
from ..store.utils import OutputSource
from ..utils.echo import echo

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
        self.logger.info("You are running the app with a standard runner.")
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
        self.input_shape = "Single"
        self.logger.debug("This is the input type: {0}".format(self.input_type))
        self.encoder = self.get_identifier_object_by_input_type()

        self.columns_info = ColumnsInformation(
            model_id=model_id, api_name=DEFAULT_API_NAME
        ).load()
        self.input_header = self.get_input_header()
        self.output_header = self.get_output_header()

        self.generic_adapter = GenericOutputAdapter(model_id, self.columns_info)

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
            if field == "Input Shape":
                self.logger.debug("Assuming input shape is Single")
                return "Single"

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
        if self.input_type[0] == "Compound":
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
        if len(self.columns_info["name"]) > 0:
            return True
        return False

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

    def get_input_header(self):
        """
        Get the input header.

        Returns
        -------
        list | None
            The input header as a list of column names, or None if the header could not be determined.
        """
        return ["key", "input"]

    def get_output_header(self):
        """
        Get the header from the columns file.

        Returns
        -------
        List | None
           Returns the header which is a list of column names, or None if the header could not be determined.
        """
        return self.columns_info.get("name", None)

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
        data_list = self.get_list_from_csv(input_data)
        data_list = [data for data in data_list]
        json_data = await self.encoder.encode_batch(data_list)
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
        data_list = []
        with open(input_data, mode="r") as file:
            reader = csv.DictReader(file)
            header = reader.fieldnames
            key = header[0] if len(header) == 1 else header[1]
            for row in reader:
                data = row.get(key)
                data_list.append(data)
        return data_list

    def serialize_to_json(self, input_data):
        """
            Serialize input data to JSON format.

            Parameters
        ----------
            input_data : str | list
                Input data which can be a file path, an input string, or a list of input strings.

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
            else:
                raise ValueError(
                    "More than one column found in input! This is not standard."
                )
        else:
            raise ValueError(
                "Input must be either a file path (string), a input string, or a list of input strings."
            )

    def is_amenable(self, output):
        """
        Check if the input and output data are amenable for a standard run.

        Parameters
        ----------
        output_data : str
            The output data to check.

        Returns
        -------
        bool
            True if the request is amenable for a standard run, False otherwise.
        """
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
        self.logger.debug("It seems amenable for standard run")
        return True

    def _serialize_output(self, result, output, model_id, api_name):
        self.logger.debug("Serializing output with generic adapter")
        self.generic_adapter._adapt_generic(result, output, model_id, api_name)
        self.logger.debug("Output serialized with generic adapter")

    def _post_batch(self, url, input_batch):
        if not input_batch:
            return []

        try:
            response = requests.post(url, json=input_batch)
            self.logger.debug("Status code: {0}".format(response.status_code))
            response.raise_for_status()
            result = response.json()
            return result
        except Exception as e:
            self.logger.error(
                "Error processing batch of size {}: {}".format(len(input_batch), e)
            )
            if len(input_batch) == 1:
                if self.output_header is not None:
                    empty_values = [None] * len(self.output_header)
                    empty_values = [
                        {k: v for k in self.output_header for v in empty_values}
                    ]
                    return empty_values
                else:
                    return [None]
            mid = len(input_batch) // 2
            left_results = self._post_batch(url, input_batch[:mid])
            right_results = self._post_batch(url, input_batch[mid:])
            return left_results + right_results

    def post(self, input, output, batch_size, output_source):
        """
        Post input data to the API in batches and get the output data.
        Finally, create a CSV file with the combined results.

        Parameters
        ----------
        input : str | list
            Input data which can be a file path, a input string, or a list of input strings.
        output : str
            Path to the output CSV file.
        batch_size : int
            Number of items per batch.
        output_source : OutputSource, optional
            Source of the output data. Default is OutputSource.LOCAL_ONLY.

        Returns
        -------
        str | None
            Path to the output CSV file if successful, None otherwise.
        """
        input_data = self.serialize_to_json(input)
        if OutputSource.is_precalculation_enabled(output_source):
            store = InferenceStoreApi(
                model_id=self.model_id, output=output, output_source=output_source
            )
            return store.get_precalculations(input_data)

        url = f"{self.url}/{self.api_name}"

        if not isinstance(input_data, list):
            input_data = [input_data]

        st = time.perf_counter()

        self.logger.info("Waiting for the server to respond...")
        results, meta = self._fetch_result(input_data, url, batch_size)
        self.logger.info("The server has responded")
        self.logger.info("Standardizing output...")
        results = self._standardize_output(input_data, results, meta)
        self.logger.debug(f"Results (chunked string): {results}"[:200])
        et = time.perf_counter()
        self.logger.info(f"All batches processed in {et - st:.4f} seconds")
        self.logger.debug("Serializing output...")
        self._serialize_output(
            json.dumps(results), output, self.model_id, self.api_name
        )
        self.logger.debug("Output serialized")
        self.logger.debug("Checking same row counts")
        matchs = self._same_row_count(input, results)
        if not matchs:
            raise Exception(
                "Inputs and outputs are not matching! Please refrain from using the results. Try again, removing bad inputs!"
            )
        ft = time.perf_counter()
        self.logger.info(f"Output is being generated within: {ft - st:.5f} seconds")
        echo(f"Output is being generated within: {ft - st:.5f} seconds")
        return output

    def _fetch_result(self, input_data, url, batch_size):
        total = len(input_data)
        overall_results = []
        meta = None

        for i in range(0, total, batch_size):
            batch_items = input_data[i : i + batch_size]
            payload = [d["input"] for d in batch_items]

            bidx = i // batch_size + 1
            echo(f"Running batch {bidx}")
            st = time.perf_counter()

            resp = self._post_batch(url, payload)

            et = time.perf_counter()
            msg = f"Batch {bidx} response fetched within: {et - st:.4f} seconds"
            self.logger.info(msg)
            echo(msg)
            try:
                if isinstance(resp, dict):
                    data = resp.json()
                else:
                    data = resp
            except ValueError:
                self.logger.error(f"Batch {bidx} returned non-JSON: {resp.text[:200]}")
                continue
            if isinstance(data, list):
                overall_results.extend(data)
            elif isinstance(data, dict):
                if "result" in data:
                    self.logger.warning("Result is in batch")
                    overall_results.extend(data["result"])
                    if "meta" in data:
                        meta = data["meta"]
                else:
                    overall_results.append(data)
            else:
                self.logger.error(f"Unexpected payload type: {type(data).__name__}")

        return overall_results, meta

    def _normalize_values(self, out):
        if isinstance(out, list):
            return out
        if out is None:
            return [None] * len(self.output_header) if self.output_header else [None]
        return list(out.values())

    def _select_key_candidates(self, out, default_keys):
        if out is None:
            return self.output_header or default_keys
        if isinstance(out, dict) and "outcome" not in out:
            return list(out.keys())
        return default_keys

    def _standardize_output(self, input_data, results, meta):
        results = list(results)
        default_keys = self._generate_default_keys(meta)

        standardized = []
        for inp, out in zip(input_data, results):
            values = self._normalize_values(out)
            key_candidates = self._select_key_candidates(out, default_keys)
            keys_flat = (
                list(key_candidates[0])
                if self._contains_dict_keys(key_candidates)
                else key_candidates
            )
            standardized.append({"input": inp, "output": dict(zip(keys_flat, values))})
        return standardized

    def _generate_default_keys(self, meta):
        if meta is None:
            return ["outcome"]

        outcomes = meta.get("outcome", [])
        if not outcomes:
            return ["outcome"]

        if len(outcomes) == 1 or "outcome" in outcomes[0]:
            return [outcomes[0]]

        return outcomes

    def _contains_dict_keys(self, lst):
        dict_keys_type = type({}.keys())
        return any(isinstance(x, dict_keys_type) for x in lst)

    def _same_row_count(self, inputs, results, include_header=False):
        def cnt(f):
            with open(f, newline="") as fp:
                n = sum(1 for _ in csv.reader(fp))
            return n if include_header else max(n - 1, 0)

        return cnt(inputs) == len(results)
