import os
import json
import collections
from pathlib import Path

from . import BaseAction
from .... import ErsiliaModel
from ....io.input import ExampleGenerator
from ....io.pure import PureDataTyper
from ....default import API_SCHEMA_FILE, MODEL_SIZE_FILE

N = 3


class ModelSniffer(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.logger.debug("Initializing model for inferring its structure")
        self.model = ErsiliaModel(
            model_id, config_json=config_json, fetch_if_not_available=False
        )
        self.logger.debug("Model successfully initialized in sniffer")
        eg = ExampleGenerator(model_id, config_json=config_json)
        self.inputs = eg.example(N, file_name=None, simple=True)
        self.logger.debug("Inputs sampled: {0}".format(len(self.inputs)))

    @staticmethod
    def __dicts_are_identical(dicts):
        for d1 in dicts:
            for d2 in dicts:
                if d1 != d2:
                    return False
        return True

    @staticmethod
    def _get_directory_size(dir):
        root_directory = Path(dir)
        bytes = sum(
            f.stat().st_size for f in root_directory.glob("**/*") if f.is_file()
        )
        return bytes

    def _get_size_in_mb(self):
        dest_dir = self._model_path(self.model_id)
        repo_dir = self._get_bundle_location(self.model_id)
        size = self._get_directory_size(dest_dir) + self._get_directory_size(repo_dir)
        mbytes = size / (1024 ** 2)
        return mbytes

    def _get_schema(self, results):
        input_schema = collections.defaultdict(list)
        output_schema = collections.defaultdict(list)
        for res in results:
            inp = res["input"]
            for k, v in inp.items():
                pdt = PureDataTyper(v)
                input_schema[k] += [pdt.get_type()]
            out = res["output"]
            if out is None:
                continue
            for k, v in out.items():
                pdt = PureDataTyper(v)
                output_schema[k] += [pdt.get_type()]
        input_schema_ = {}
        for k, v in input_schema.items():
            if self.__dicts_are_identical(v):
                input_schema_[k] = v[0]
            else:
                self.looger.error("Input data types are not consistent")
        output_schema_ = {}
        for k, v in output_schema.items():
            if self.__dicts_are_identical(v):
                output_schema_[k] = v[0]
            else:
                self.logger.error("Output data types are not consistent")
        meta = self.model.autoservice._latest_meta
        for k, v in output_schema_.items():
            output_schema_[k]["meta"] = meta[k]
        schema = {"input": input_schema_, "output": output_schema_}
        return schema

    def sniff(self):
        self.logger.debug("Sniffing model")
        self.logger.debug("Getting model size")
        size = self._get_size_in_mb()
        self.logger.debug("Mode size is {0} MB".format(size))
        path = os.path.join(self._model_path(self.model_id), MODEL_SIZE_FILE)
        with open(path, "w") as f:
            json.dump({"size": size, "units": "MB"}, f, indent=4)
        self.logger.debug("Serving model")
        self.model.autoservice.serve()
        self.logger.debug("Iterating over APIs")
        all_schemas = {}
        for api_name in self.model.autoservice.get_apis():
            self.logger.debug("Running API: {0}".format(api_name))
            self.logger.debug(self.inputs)
            results = [
                result for result in self.model.autoservice.api(api_name, self.inputs)
            ]
            schema = self._get_schema(results)
            self.logger.debug(schema)
            all_schemas[api_name] = schema
        path = os.path.join(self._model_path(self.model_id), API_SCHEMA_FILE)
        with open(path, "w") as f:
            json.dump(all_schemas, f, indent=4)
        self.logger.debug("API schema saved at {0}".format(path))
