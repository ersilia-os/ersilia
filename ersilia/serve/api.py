import os
import csv
import requests
import json
import collections
import tempfile
import time

from ..io.input import GenericInputAdapter
from ..io.output import GenericOutputAdapter
from ..lake.interface import IsauraInterface
from .. import logger
from .. import ErsiliaBase
from .schema import ApiSchema

from ..utils.exceptions_utils.api_exceptions import InputFileNotFoundError


class Api(object):
    def __init__(self, model_id, url, api_name, save_to_lake, config_json):
        self.config_json = config_json
        self.model_id = model_id
        self.input_adapter = GenericInputAdapter(self.model_id, config_json=config_json)
        self.output_adapter = GenericOutputAdapter(config_json=config_json)
        self.lake = IsauraInterface(
            model_id=model_id, api_name=api_name, config_json=config_json
        )
        self.save_to_lake = save_to_lake
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

    def _do_post(self, input, output):
        url = "{0}/{1}".format(self.url, self.api_name)
        if self._do_sleep:
            time.sleep(3)
        response = requests.post(url, json=input)
        self.logger.debug("Status code: {0}".format(response.status_code))
        if response.status_code == 200:
            result_ = response.json()
            result_ = self.output_adapter.refactor_response(result_)
            result = []
            for i, o in zip(input, result_):
                result += [{"input": i, "output": o}]
            result = json.dumps(result, indent=4)
            return self.__result_returner(result, output)
        else:
            self.logger.error("Status Code: {0}".format(response.status_code))
            return None

    def _post(self, input, output):
        result = self._do_post(input, output)
        if result is None and self._batch_size > 1 and len(input) > 1:
            # if batch predictions didn't work, do one by one
            self.logger.warning(
                "Batch prediction didn't seem to work. Doing predictions one by one..."
            )
            result = []
            for inp_one in input:
                r = self._do_post([inp_one], output=None)
                if r is None:
                    r = [{"input": inp_one, "output": self._empty_output}]
                else:
                    r = json.loads(r)
                result += r
            result = json.dumps(result, indent=4)
            result = self.__result_returner(result, output)
        if result is None and len(input) == 1:
            # if the only one prediction did not work, return empty
            result = [{"input": input[0], "output": self._empty_output}]
            result = json.dumps(result, indent=4)
            result = self.__result_returner(result, output)
        return result

    def post_only_calculations(self, input, output, batch_size):
        self._batch_size = batch_size
        if output is not None:
            tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
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
            for input in self.input_adapter.adapt(input, batch_size=batch_size):
                result = json.loads(self._post(input, output))
                for r in result:
                    yield r

    def _post_reads(self, input, output):
        results = self.lake.read(input)
        results = json.dumps(results, indent=4)
        return self.__result_returner(results, output)

    def post_only_reads(self, input, output, batch_size):
        self._batch_size = batch_size
        if output is not None:
            tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
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

    def post_amenable_to_h5(self, input, output, batch_size):
        self.logger.debug(
            "Checking for already available calculations in the data lake"
        )
        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        done_input = os.path.join(tmp_folder, "done_input.csv")
        todo_input = os.path.join(tmp_folder, "todo_input.csv")
        cur_idx = 0
        for input in self.input_adapter.adapt(input, batch_size=batch_size):
            self.logger.debug("Inspecting {0}...".format(cur_idx))
            done, todo = self.lake.done_todo(input)
            self._write_done_todo_file(cur_idx, done_input, done)
            self._write_done_todo_file(cur_idx, todo_input, todo)
            cur_idx += len(done) + len(todo)
        if self.__is_empty_file(done_input):
            done_output = None
        else:
            done_output = os.path.join(tmp_folder, "done_output.json")
        if self.__is_empty_file(todo_input):
            todo_output = None
        else:
            todo_output = os.path.join(tmp_folder, "todo_output.json")
        if done_output is not None:
            self.logger.debug("Reading from data well of {0}".format(self.model_id))
            for _ in self.post_only_reads(
                input=done_input, output=done_output, batch_size=batch_size
            ):
                continue
        if todo_output is not None:
            self.logger.debug("Calculating using model {0}".format(self.model_id))
            self.logger.debug("Saving in {0}".format(todo_output))
            for _ in self.post_only_calculations(
                input=todo_input, output=todo_output, batch_size=batch_size
            ):
                continue

            with open(todo_output, "r") as f:
                results = json.load(f)
            if self.save_to_lake:
                self.logger.debug("Saving calculations in the lake")
                self.lake.write(results)

        self.logger.debug("Rearranging and returning")
        results = []
        for result in self._process_done_todo_results(
            done_input, todo_input, done_output, todo_output
        ):
            results += [result]
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

    def _unique_input(self, input):
        mapping = collections.defaultdict(list)
        unique_input = []
        for i, inp in enumerate(self.input_adapter.adapt_one_by_one(input)):
            key = inp["key"]
            if key not in mapping:
                unique_input += [inp]
            mapping[key] += [i]
        return unique_input, mapping

    def post_unique_input(self, input, output, batch_size):
        schema = ApiSchema(model_id=self.model_id, config_json=self.config_json)
        if (
            not schema.isfile()
            or not schema.is_h5_serializable(api_name=self.api_name)
            or not self.lake.is_available
        ):
            for res in self.post_only_calculations(input, output, batch_size):
                yield res
        else:
            for res in self.post_amenable_to_h5(input, output, batch_size):
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

    def post(self, input, output, batch_size):
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
        return self.output_adapter.meta()
