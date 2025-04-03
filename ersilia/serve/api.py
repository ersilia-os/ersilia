import collections
import csv
import json
import os
import time

import requests

from .. import ErsiliaBase, logger
from ..default import API_SCHEMA_FILE, DEFAULT_API_NAME
from ..io.input import GenericInputAdapter
from ..io.output import GenericOutputAdapter
from ..utils.exceptions_utils.api_exceptions import InputFileNotFoundError
from ..utils.logging import make_temp_dir
from .schema import ApiSchema


class Api(ErsiliaBase):
    """
    Class to interact with the API for a given model.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    url : str
        The URL of the API.
    api_name : str
        The name of the API.
    config_json : dict
        Configuration in JSON format.

    Examples
    --------
    .. code-block:: python

        api = Api(
            model_id="eosxxxx",
            url="http://0.0.0.0:25512/",
            api_name="run",
            config_json={},
        )
        result = api.post(
            input="input.json",
            output="output.csv",
            batch_size=10,
        )
    """

    def __init__(self, model_id, url, api_name, config_json):
        ErsiliaBase.__init__(self, config_json=None)
        self.config_json = config_json
        self.model_id = model_id
        self.input_adapter = GenericInputAdapter(
            model_id=self.model_id, config_json=config_json
        )
        self.output_adapter = GenericOutputAdapter(
            model_id=self.model_id, config_json=config_json
        )
        if url[-1] == "/":
            self.url = url[:-1]
        else:
            self.url = url
        self.api_name = api_name
        self.logger = logger
        self.logger.debug(
            "API {0}:{1} initialized at URL {2}".format(model_id, api_name, url)
        )
        try:
            self._empty_output = ApiSchema(
                model_id=model_id, config_json=None
            ).empty_output_by_api(api_name=api_name)
        except:
            self.logger.info("No empty output available")
            self._empty_output = None
        if self._is_during_fetch():
            self._do_sleep = True
        else:
            self._do_sleep = False

    def _is_during_fetch(self):
        base = ErsiliaBase(config_json=self.config_json, credentials_json=None)
        path = os.path.join(
            base._get_bundle_location(model_id=self.model_id), "status.json"
        )
        if not os.path.exists(path):
            return True
        with open(path, "r") as f:
            data = json.load(f)
        if data["done"]:
            return False
        else:
            return True

    def __result_returner(self, result, output):
        if output is None:
            return self.output_adapter.adapt(
                result, output, model_id=self.model_id, api_name=self.api_name
            )
        else:
            self.logger.debug("Working on output: {0}".format(output))
            self.output_adapter.adapt(
                result, output, model_id=self.model_id, api_name=self.api_name
            )
            return [{"output": output}]

    def _post_batch(self, url, input_batch):
        if not input_batch:
            return []

        try:
            response = requests.post(url, json=input_batch)

            self.logger.debug("Status code: {0}".format(response.status_code))
            response.raise_for_status()
            result = response.json()
            result = self.output_adapter.refactor_response(result)
            return result
        except Exception as e:
            self.logger.error(
                "Error processing batch of size {}: {}".format(len(input_batch), e)
            )
            if len(input_batch) == 1:
                return [self._empty_output]
            mid = len(input_batch) // 2
            left_results = self._post_batch(url, input_batch[:mid])
            right_results = self._post_batch(url, input_batch[mid:])
            return left_results + right_results

    def _do_post(self, input, output):
        url = "{0}/{1}".format(self.url, self.api_name)
        if self._do_sleep:
            time.sleep(3)

        results = self._post_batch(url, input)

        combined_results = []
        for compound, res in zip(input, results):
            combined_results.append({"input": compound, "output": res})

        result_json = json.dumps(combined_results, indent=4)
        return self.__result_returner(result_json, output)

    def _post(self, input, output):
        import time

        st = time.time()
        result = self._do_post(input, output)
        result = self._standardize_schema(result)
        et = time.time()
        self.logger.debug("Time taken: {0}".format(et - st))
        return result

    def post(self, input, output, batch_size):
        """
        Post input data to the API and get the result.

        Parameters
        ----------
        input : str
            The input data file or data.
        output : str
            The output data file.
        batch_size : int
            The batch size for processing.

        Yields
        ------
        dict
            The result of the API call.
        """
        if self._is_input_file(input):
            if not os.path.exists(input):
                raise InputFileNotFoundError(file_name=input)
        self.logger.debug("Posting to {0}".format(self.api_name))
        self.logger.debug("Batch size {0}".format(batch_size))
        unique_input, mapping = self._unique_input(input)
        results_ = {}
        for res in self.post_unique_input(
            input=unique_input, output=None, batch_size=batch_size
        ):
            for i in mapping[res["input"]["key"]]:
                results_[i] = res
        self.logger.debug("Done with unique posting")
        sorted_idxs = sorted(results_.keys())
        results = [results_[i] for i in sorted_idxs]
        if output is not None:
            results = json.dumps(results)
            self.output_adapter.adapt(
                results, output, model_id=self.model_id, api_name=self.api_name
            )
            for o in [output]:
                yield o
        else:
            for result in results:
                yield result

    def meta(self):
        """
        Get metadata from the output adapter.

        Returns
        -------
        dict
            Metadata information.
        """
        return self.output_adapter.meta()

    def post_only_calculations(self, input, output, batch_size):
        """
        Post input data to the API and get the result, performing only calculations.

        Parameters
        ----------
        input : str
            The input data file or data.
        output : str
            The output data file.
        batch_size : int
            The batch size for processing.

        Yields
        ------
        dict
            The result of the API call.
        """
        self._batch_size = batch_size
        if output is not None:
            tmp_folder = make_temp_dir(prefix="ersilia-")
            fmt = output.split(".")[-1]
            output_base = ".".join(os.path.basename(output).split(".")[:-1])
            i = 0
            subfiles = []
            for input in self.input_adapter.adapt(input, batch_size=batch_size):
                subfile = os.path.join(
                    tmp_folder, "{0}-chunk-{1}.{2}".format(output_base, i, fmt)
                )
                self._post(input, subfile)
                subfiles += [subfile]
                i += 1
            self.output_adapter.merge(subfiles, output)
            for o in [output]:
                yield o
        else:
            self.logger.debug("Posting only calculations")
            for input in self.input_adapter.adapt(input, batch_size=batch_size):
                result = json.loads(self._post(input, output))
                for r in result:
                    yield r

    def post_only_reads(self, input, output, batch_size):
        """
        Post input data to the API and get the result, performing only reads.

        Parameters
        ----------
        input : str
            The input data file or data.
        output : str
            The output data file.
        batch_size : int
            The batch size for processing.

        Yields
        ------
        dict
            The result of the API call.
        """
        self._batch_size = batch_size
        if output is not None:
            tmp_folder = make_temp_dir(prefix="ersilia-")
            fmt = output.split(".")[-1]
            output_base = ".".join(os.path.basename(output).split(".")[:-1])
            i = 0
            subfiles = []
            for input in self.input_adapter.adapt(input, batch_size=batch_size):
                subfile = os.path.join(
                    tmp_folder, "{0}-chunk-{1}.{2}".format(output_base, i, fmt)
                )
                self._post_reads(input, subfile)
                subfiles += [subfile]
                i += 1
            self.output_adapter.merge(subfiles, output)
            for o in [output]:
                yield o
        else:
            for input in self.input_adapter.adapt(input, batch_size=batch_size):
                result = json.loads(self._post_reads(input, output))
                for r in result:
                    yield r

    def post_unique_input(self, input, output, batch_size):
        """
        Post unique input data to the API and get the result.

        Parameters
        ----------
        input : str
            The input data file or data.
        output : str
            The output data file.
        batch_size : int
            The batch size for processing.

        Yields
        ------
        dict
            The result of the API call.
        """
        schema = ApiSchema(model_id=self.model_id, config_json=self.config_json)
        self.logger.debug(" Post unique input data to the API")
        if schema.isfile():
            for res in self.post_only_calculations(input, output, batch_size):
                yield res

    def _is_input_file(self, input):
        if type(input) is str:
            if input.endswith(".csv"):
                return True
            if input.endswith(".tst"):
                return True
            if input.endswith(".json"):
                return True
            if input.endswith(".txt"):
                return True
        return False

    def _unique_input(self, input):
        mapping = collections.defaultdict(list)
        unique_input = []
        for i, inp in enumerate(self.input_adapter.adapt_one_by_one(input)):
            key = inp["key"]
            if key not in mapping:
                unique_input += [inp]
            mapping[key] += [i]
        return unique_input, mapping

    def _write_done_todo_file(self, cur_idx, filename, data):
        with open(filename, "a+") as f:
            writer = csv.writer(f)
            for d in data:
                idx = cur_idx + d["idx"]
                key = d["key"]
                inp = d["input"]
                txt = d["text"]
                writer.writerow([idx, key, inp, txt])

    def _process_done_todo_results(
        self, done_input, todo_input, done_output, todo_output
    ):
        mapping = {}
        if done_output is not None:
            with open(done_input, "r") as f:
                reader = csv.reader(f)
                for i, r in enumerate(reader):
                    mapping[int(r[0])] = (i, True)
            with open(done_output, "r") as f:
                done_output_data = json.load(f)
        else:
            done_output_data = {}
        if todo_output is not None:
            with open(todo_input, "r") as f:
                reader = csv.reader(f)
                for i, r in enumerate(reader):
                    mapping[int(r[0])] = (i, False)
            with open(todo_output, "r") as f:
                todo_output_data = json.load(f)
        else:
            todo_output_data = {}
        for j in range(len(mapping)):
            i, is_done = mapping[j]
            if is_done:
                yield done_output_data[i]
            else:
                yield todo_output_data[i]

    @staticmethod
    def __is_empty_file(filename):
        if os.stat(filename).st_size == 0:
            return True
        else:
            return False

    def _load_api_schema(self):
        outcome_key, reference_keys = "outcome", []
        self.schema_file = os.path.join(
            self._model_path(self.model_id), API_SCHEMA_FILE
        )
        if os.path.exists(self.schema_file):
            with open(self.schema_file, "r") as f:
                schema = json.load(f)

            output_schema = schema.get(DEFAULT_API_NAME, {}).get("output", {})
            if not output_schema:
                raise ValueError("No output schema found in the API schema file.")

            outcome_key, outcome_details = next(iter(output_schema.items()))
            reference_keys = outcome_details.get("meta", [])

        return outcome_key, reference_keys

    def _detect_mismatch(self, data):
        outcome_key, reference_keys = self._load_api_schema()
        mismatches = []
        expected_keys = set(reference_keys)

        for i, entry in enumerate(data):
            output = entry.get("output", {})

            if outcome_key in output and isinstance(output[outcome_key], dict):
                container = output[outcome_key]
            else:
                container = output

            if isinstance(container, list):
                container = {
                    key: container[idx] if idx < len(container) else None
                    for idx, key in enumerate(reference_keys)
                }

            if isinstance(container, dict):
                actual_keys = set(container.keys())
            else:
                actual_keys = set()

            missing_keys = list(expected_keys - actual_keys)
            extra_keys = list(actual_keys - expected_keys)

            if missing_keys or extra_keys:
                mismatches.append(
                    {
                        "index": i,
                        "missing_keys": missing_keys,
                        "extra_keys": extra_keys,
                    }
                )

        return mismatches

    def _standardize_schema(self, data, force_list=True):
        data = json.loads(data)
        outcome_key, reference_keys = self._load_api_schema()

        if force_list:
            for entry in data:
                output = entry.get("output", {})
                container = self._extract_data_container_force(
                    output, outcome_key, reference_keys
                )
                standardized = {key: container.get(key, None) for key in reference_keys}
                entry["output"] = {
                    outcome_key: [standardized[key] for key in reference_keys]
                }
        else:
            mismatches = self._detect_mismatch(data)
            for mismatch in mismatches:
                idx = mismatch["index"]
                entry = data[idx]
                output = entry.get("output", {})

                if outcome_key in output and isinstance(output[outcome_key], dict):
                    container = output[outcome_key]
                    container_is_nested = True
                else:
                    container = output
                    container_is_nested = False

                container = self._normalize_container(container, reference_keys)
                if container_is_nested:
                    output[outcome_key] = [
                        container.get(key, None) for key in reference_keys
                    ]
                else:
                    entry["output"] = container

        return json.dumps(data, indent=4)

    def _extract_data_container_force(self, output, outcome_key, reference_keys):
        if outcome_key in output:
            nested = output[outcome_key]
            if isinstance(nested, dict):
                return nested
            elif isinstance(nested, list):
                return self._convert_list_to_dict(nested, reference_keys)
            else:
                return {reference_keys[0]: nested} if reference_keys else {}
        else:
            if isinstance(output, list):
                return self._convert_list_to_dict(output, reference_keys)
            elif isinstance(output, dict):
                return output
            else:
                return {reference_keys[0]: output} if reference_keys else {}

    def _normalize_container(self, container, reference_keys):
        if isinstance(container, list):
            return self._convert_list_to_dict(container, reference_keys)
        elif isinstance(container, dict):
            if not any(key in reference_keys for key in container.keys()):
                values = list(container.values())
                return {
                    key: values[i] if i < len(values) else None
                    for i, key in enumerate(reference_keys)
                }
            else:
                new_container = {}
                extra_values = [
                    container[k] for k in container if k not in reference_keys
                ]
                for key in reference_keys:
                    if key in container:
                        new_container[key] = container[key]
                    elif extra_values:
                        new_container[key] = extra_values.pop(0)
                    else:
                        new_container[key] = None
                return new_container
        else:
            return {reference_keys[0]: container} if reference_keys else container

    def _convert_list_to_dict(self, lst, reference_keys):
        return {
            key: lst[i] if i < len(lst) else None
            for i, key in enumerate(reference_keys)
        }
