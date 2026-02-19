import asyncio
import csv
import importlib
import json
import os
import time
from collections import Counter

import nest_asyncio
import pandas as pd
import requests

from .. import ErsiliaBase
from ..core.session import Session
from ..default import (
    DEFAULT_API_NAME,
    EXAMPLE_STANDARD_INPUT_CSV_FILENAME,
    EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME,
    INFORMATION_FILE,
)
from ..hub.content.columns_information import ColumnsInformation
from ..io.output import GenericOutputAdapter
from ..store.isaura import IsauraStore
from ..utils.echo import echo, spinner
from ..utils.ports import _ensure_ready, normalize_connect_url

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
        self.model_id = model_id
        if url[-1] == "/":
            self.url = url[:-1]
        else:
            self.url = url
        self.url = normalize_connect_url(self.url)
        self.api_name = DEFAULT_API_NAME
        self.path = os.path.abspath(self._model_path(self.model_id))
        metadata = self._read_information_file()
        self.standard_input_csv = os.path.join(
            self.path, EXAMPLE_STANDARD_INPUT_CSV_FILENAME
        )
        self.standard_output_csv = os.path.join(
            self.path, EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME
        )
        self.input_type = self._read_field_from_metadata(metadata, "Input")
        self.input_shape = "Single"
        self.encoder = self.get_identifier_object_by_input_type()
        self.columns_info = ColumnsInformation(
            model_id=model_id, api_name=DEFAULT_API_NAME
        ).load()
        self.input_header = self.get_input_header()
        self.output_header = self.get_output_header()
        self.generic_adapter = GenericOutputAdapter(model_id, self.columns_info)
        self.isaura_store = IsauraStore()
        self.session = Session(config_json=config_json)
        store_info = self.session.current_store_status()
        self.local_cache = store_info[-1]
        self.write_store = store_info[1]
        self.read_store = store_info[0]
        echo("Standard API runner initialized", fg="green")

    def _read_information_file(self):
        try:
            return json.load(open(os.path.join(self.path, INFORMATION_FILE), "r"))
        except:
            return None

    def _read_field_from_metadata(self, meta, field):
        if not meta:
            return None
        if "metadata" in meta and field in meta["metadata"]:
            return meta["metadata"][field]
        if "card" in meta and field in meta["card"]:
            return meta["card"][field]
        return None

    def get_identifier_object_by_input_type(self):
        """
        Get the identifier object by input type.

        Returns
        -------
        object
            The identifier object.
        """
        module_path = f"ersilia.utils.identifiers.{self.input_type[0].lower()}"
        return importlib.import_module(module_path).Identifier()

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
        return self.input_type[0] == "Compound"

    def is_output_type_standardizable(self):
        """
        Check if the output type is standardizable.

        Returns
        -------
        bool
            True if the output type is standardizable, False otherwise.
        """
        return len(self.columns_info["name"]) > 0

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
        return isinstance(output_data, str) and output_data.endswith(".csv")

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
            The header which is a list of column names, or None if the header could not be determined.
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
        echo("Serializing CSV input", fg="cyan")
        json_data = []
        with open(input_data, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                key = self.encoder.encode(row[0])
                json_data.append({"key": key, "input": row[0], "text": row[0]})
        echo("CSV serialization complete", fg="green")
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
        data_list = [d for d in data_list]
        result = await self.encoder.encode_batch(data_list)
        return result

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
                data_list.append(row.get(key))
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
            return False
        if not self.is_output_type_standardizable():
            return False
        return True

    def _serialize_output(self, result, output, model_id, api_name):
        _, df = self.generic_adapter._adapt_generic(result, output, model_id, api_name)
        echo("Output serialization complete", fg="green")
        self.logger.info("Output serialization complete")
        return df

    def _missing_inputs(self, check_dict, inputs):
        return [m for m in inputs if not check_dict[m]]

    def _found_inputs(self, check_dict, inputs):
        return [m for m in inputs if check_dict[m]]

    def _post_batch(self, url, input_batch):
        if not input_batch:
            return []

        params = {
            "fetch_cache": self.local_cache,
            "save_cache": self.local_cache,
        }

        def do_request(batch):
            response = requests.post(url, params=params, json=batch)
            response.raise_for_status()
            return response.json()

        return spinner(
            f"Fetching batch of size {len(input_batch)}", do_request, input_batch
        )

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
        echo("Preparing input data", fg="cyan")
        input_data = self.serialize_to_json(input)
        url = f"{self.url}/{self.api_name}"

        if not isinstance(input_data, list):
            input_data = [input_data]

        st = time.perf_counter()
        _ensure_ready(self=self, root=self.url)
        echo("Waiting for server response", fg="cyan")
        self.logger.debug("Waiting for server response")
        results, meta = self._fetch_result(input_data, url, batch_size)
        echo("Server response received", fg="green")
        self.logger.info("Server response received")

        results = self._standardize_output(input_data, results, meta)
        et = time.perf_counter()
        echo(f"All batches processed in {et - st:.4f} seconds", fg="green")
        self.logger.info(f"All batches processed in {et - st:.4f} seconds")

        echo("Finalizing the output", fg="cyan")
        self.logger.debug("Finalizing the output")
        df = self._serialize_output(
            json.dumps(results), output, self.model_id, self.api_name
        )
        if self.write_store:
            df = pd.DataFrame(
                data=df.data, columns=["key", "input"] + self.output_header, dtype=str
            )
            echo("Writing results to Isaura store", fg="cyan")
            self.logger.info("Writing results to Isaura store")
            self.isaura_store.write(df=df)

        matchs = self._same_row_count(input_data, results)
        if not matchs:
            raise Exception("Inputs and outputs are not matching")

        ft = time.perf_counter()
        echo(f"Output generated in {ft - st:.5f} seconds", fg="green")
        self.logger.info(f"Output generated to {output} in {ft - st:.5f} seconds")

        return output

    def _check_found_and_key(self, check_dict, smi):
        if check_dict is None:
            return False, None

        if isinstance(check_dict, dict):
            v = check_dict.get(smi)
            if v is None:
                return False, None
            if isinstance(v, bool):
                return v, None
            if isinstance(v, str):
                return True, v
            if isinstance(v, dict):
                if "found" in v:
                    return bool(v.get("found")), v.get("key")
                if "exists" in v:
                    return bool(v.get("exists")), v.get("key")
                if "key" in v:
                    return True, v.get("key")
                return bool(v), v.get("key")
            if isinstance(v, (list, tuple)):
                if len(v) == 0:
                    return False, None
                if isinstance(v[0], bool):
                    k = v[1] if len(v) > 1 and isinstance(v[1], str) else None
                    return bool(v[0]), k
                if isinstance(v[0], str):
                    return True, v[0]
                return True, None
            return True, None

        if isinstance(check_dict, (set, list, tuple)):
            return smi in check_dict, None

        return False, None

    def _fetch_result(self, input_data, url, batch_size):
        total = len(input_data)
        meta = None
        overall_results = [None] * total
        all_replacements = []

        def write_replacements_txt(repls):
            out_path = getattr(self, "output_path", None) or "output.csv"
            base, _ = os.path.splitext(out_path)
            path = base + ".ann_replacements.txt"

            counts = Counter(
                (r["input"], r["lookup_input"], r.get("source", "")) for r in repls
            )

            with open(path, "w", encoding="utf-8") as f:
                f.write("count\tinput\tlookup_input\tsource\n")
                for (inp, lk, src), c in counts.most_common():
                    f.write(f"{c}\t{inp}\t{lk}\t{src}\n")

            return path

        for i in range(0, total, batch_size):
            batch_items = input_data[i : i + batch_size]

            payload = [d["input"] for d in batch_items]
            lookup_payload = [d.get("lookup_input") or d["input"] for d in batch_items]

            if IsauraStore.is_installed() and self.read_store and lookup_payload:
                check_dict = self.isaura_store.check(lookup_payload)
                found = []
                missed = []
                for smi in lookup_payload:
                    is_found, _ = self._check_found_and_key(check_dict, smi)
                    (found if is_found else missed).append(smi)
            else:
                found = []
                missed = lookup_payload

            found_u = list(dict.fromkeys(found))
            missed_u = list(dict.fromkeys(missed))

            if found_u and IsauraStore.is_installed():
                cache_df = self.isaura_store.read(found_u)
            else:
                cache_df = pd.DataFrame(columns=["key", "input"])

            api_values = []
            if missed_u:
                bidx = i // batch_size + 1
                echo(f"Running batch {bidx}", fg="cyan")
                self.logger.debug(f"Running batch {bidx}")
                st = time.perf_counter()

                resp = self._post_batch(url, missed_u)

                et = time.perf_counter()
                echo(f"Batch {bidx} fetched in {et - st:.4f} seconds", fg="green")
                self.logger.info(f"Batch {bidx} fetched in {et - st:.4f} seconds")

                try:
                    data = resp.json() if hasattr(resp, "json") else resp
                except Exception:
                    data = None

                if isinstance(data, list):
                    api_values = data
                elif isinstance(data, dict) and "result" in data:
                    api_values = data["result"]
                    if isinstance(data.get("meta"), dict):
                        if meta is None:
                            meta = {}
                        meta.update(data["meta"])
                elif data is not None:
                    api_values = [data] * len(missed_u)
                else:
                    api_values = []

            batch_results, batch_repls = self._merge_cache_and_api(
                payload=payload,
                lookup_payload=lookup_payload,
                found=found_u,
                cache_df=cache_df,
                missed=missed_u,
                api_values=api_values,
            )

            all_replacements.extend(batch_repls)
            overall_results[i : i + len(batch_results)] = batch_results

        if all_replacements:
            path = write_replacements_txt(all_replacements)
            self.logger.warning(
                f"ANN replacements detected count={len(all_replacements)} file={path}"
            )
            if meta is None:
                meta = {}
            meta["ann_replacements_count"] = len(all_replacements)
            meta["ann_replacements_file"] = path

        return overall_results, meta

    def _write_ann_replacements_txt(self, replacements):
        out_path = getattr(self, "output_path", None) or "output.csv"
        base, _ = os.path.splitext(out_path)
        path = base + ".ann_replacements.txt"

        counts = Counter(
            (r["input"], r["lookup_input"], r["source"]) for r in replacements
        )
        with open(path, "w", encoding="utf-8") as f:
            f.write("count\tinput\tlookup_input\tsource\n")
            for (inp, lk, src), c in counts.most_common():
                f.write(f"{c}\t{inp}\t{lk}\t{src}\n")

        return path

    def _merge_cache_and_api(
        self, payload, lookup_payload, found, cache_df, missed, api_values
    ):
        t0 = time.perf_counter()
        value_cols = (
            [c for c in cache_df.columns if c not in ("key", "input")]
            if cache_df is not None
            else []
        )

        cache_map = {}
        if cache_df is not None and not cache_df.empty and found:
            if len(cache_df) == len(found):
                rows = (
                    cache_df[value_cols].to_dict(orient="records")
                    if value_cols
                    else [{}] * len(cache_df)
                )
                for smi, vals in zip(found, rows):
                    cache_map[smi] = vals
            else:
                self.logger.warning(
                    f"cache alignment mismatch found={len(found)} cache_rows={len(cache_df)}"
                )

        api_map = {}
        missed = missed or []
        api_values = api_values or []
        if missed:
            if not isinstance(api_values, (list, tuple)):
                api_values = [api_values] * len(missed)
            n = min(len(missed), len(api_values))
            for smi, vals in zip(missed[:n], api_values[:n]):
                if isinstance(vals, dict):
                    api_map[smi] = vals
                elif isinstance(vals, (list, tuple)):
                    api_map[smi] = {f"value_{j}": v for j, v in enumerate(vals)}
                else:
                    api_map[smi] = {"value": vals}

        results = []
        replacements = []

        for orig_smi, lookup_smi in zip(payload, lookup_payload):
            lookup_smi = (
                lookup_smi.strip() if isinstance(lookup_smi, str) else lookup_smi
            )

            if lookup_smi in cache_map:
                results.append({"input": orig_smi, **cache_map[lookup_smi]})
                if lookup_smi != orig_smi:
                    replacements.append(
                        {
                            "input": orig_smi,
                            "lookup_input": lookup_smi,
                            "source": "cache",
                        }
                    )
                continue

            if lookup_smi in api_map:
                results.append({"input": orig_smi, **api_map[lookup_smi]})
                if lookup_smi != orig_smi:
                    replacements.append(
                        {"input": orig_smi, "lookup_input": lookup_smi, "source": "api"}
                    )
                continue

            results.append({"input": orig_smi, **{v: None for v in self.output_header}})

        self.logger.debug(
            f"_merge_cache_and_api done payload={len(payload)} found={len(found or [])} missed={len(missed)} "
            f"cache_rows={len(cache_df) if cache_df is not None else 0} replacements={len(replacements)} "
            f"dt_total={(time.perf_counter() - t0):.6f}s"
        )
        return results, replacements

    def _normalize_values(self, out):
        if isinstance(out, list):
            return out
        if out is None:
            return [None] * len(self.output_header) if self.output_header else [None]
        return list(out.values())

    def _standardize_output(self, input_data, results, meta):
        results = list(results)
        standardized = []
        for inp, out in zip(input_data, results):
            values = self._normalize_values(out)
            keys_flat = self.input_header[1:] + self.output_header
            standardized.append({"input": inp, "output": dict(zip(keys_flat, values))})
        return standardized

    def _same_row_count(self, inputs, results):
        return len(inputs) == len(results)
