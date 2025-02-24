import collections
import csv
import json
import os
from pathlib import Path

from .... import ErsiliaBase, ErsiliaModel, throw_ersilia_exception
from ....default import (
    API_SCHEMA_FILE,
    MODEL_SIZE_FILE,
    PREDEFINED_EXAMPLE_FILES,
)
from ....io.annotated import AnnotatedDataTyper
from ....io.input import ExampleGenerator
from ....io.pure import PureDataTyper
from ....utils.exceptions_utils.exceptions import EmptyOutputError
from ....utils.exceptions_utils.fetch_exceptions import (
    OutputDataTypesNotConsistentError,
)
from ....utils.paths import get_metadata_from_base_dir
from . import BaseAction


class BuiltinExampleReader(ErsiliaBase):
    """
    Reads built-in examples for a BentoML model.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    config_json : dict
        Configuration settings for the reader.

    Methods
    -------
    has_builtin_example() -> bool
        Checks if a built-in example exists.
    example(n=3) -> list
        Returns a list of examples.
    """

    def __init__(self, model_id: str, config_json: dict):
        super().__init__(config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.example_file = None
        for pf in PREDEFINED_EXAMPLE_FILES:
            example_file = os.path.join(
                self._model_path(self.model_id),
                pf,
            )
            if os.path.exists(example_file):
                self.example_file = example_file
                break

    def has_builtin_example(self) -> bool:
        """
        Checks if a built-in example exists.

        Returns
        -------
        bool
            True if a built-in example exists, False otherwise.
        """
        if self.example_file is None:
            return False
        if os.path.exists(self.example_file):
            return True
        else:
            return False

    def example(self, n: int = 3) -> list:
        """
        Returns a list of examples.

        Parameters
        ----------
        n : int, optional
            Number of examples to return, by default 3.

        Returns
        -------
        list
            List of examples.
        """
        data = []
        with open(self.example_file, "r") as f:
            reader = csv.reader(f)
            next(reader)
            for r in reader:
                data += [r[0]]
        return data[:n]


class ModelSniffer(BaseAction):
    """
    Infers the structure of a model by sniffing its inputs and outputs.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    config_json : dict
        Configuration settings for the sniffer.

    Methods
    -------
    sniff()
        Infers the structure of the model.
    """

    def __init__(self, model_id: str, config_json: dict):
        super().__init__(
            model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.logger.debug("Initializing model for inferring its structure")
        self.model = ErsiliaModel(
            model_id, config_json=config_json, fetch_if_not_available=False
        )
        self.model_id = model_id
        self.logger.debug("Model successfully initialized in sniffer")
        er = BuiltinExampleReader(model_id, config_json=config_json)
        if er.has_builtin_example():
            self.logger.debug("Built-in example found")
            self.inputs = er.example()
        else:
            self.logger.debug("No built-in example available. Generating a test one.")
            eg = ExampleGenerator(model_id, config_json=config_json)
            self.inputs = eg.test()
        self.logger.debug("Inputs sampled: {0}".format(len(self.inputs)))

    @staticmethod
    def __dicts_are_identical(dicts: list) -> bool:
        for d1 in dicts:
            for d2 in dicts:
                if d1 != d2:
                    return False
        return True

    @staticmethod
    def _get_directory_size(dir: str) -> int:
        root_directory = Path(dir)
        bytes = sum(
            f.stat().st_size for f in root_directory.glob("**/*") if f.is_file()
        )
        return bytes

    def _get_size_in_mb(self) -> float:
        dest_dir = self._model_path(self.model_id)
        repo_dir = self._get_bundle_location(self.model_id)
        size = self._get_directory_size(dest_dir) + self._get_directory_size(repo_dir)
        mbytes = size / (1024**2)
        return mbytes

    def _get_output_ann_type(self):
        dest_dir = self._model_path(self.model_id)
        try:
            md = get_metadata_from_base_dir(dest_dir)
        except FileNotFoundError:
            self.logger.debug("Metadata file not available")
            return
        return md["Output Type"][0]  # TODO Account for mixed types

    def _get_output_ann_shape(self):
        dest_dir = self._model_path(self.model_id)
        try:
            md = get_metadata_from_base_dir(dest_dir)
        except FileNotFoundError:
            self.logger.debug("Metadata file not available")
            return
        return md["Output Shape"]  # TODO Account for mixed types

    @throw_ersilia_exception()
    def _get_schema(self, results: list) -> dict:
        input_schema = collections.defaultdict(list)
        output_schema = collections.defaultdict(list)
        for res in results:
            inp = res["input"]
            for k, v in inp.items():
                dt = PureDataTyper(v)
                input_schema[k] += [dt.get_type()]
            out = res["output"]
            if out is None:
                continue
            for k, v in out.items():
                ann_type = self._get_output_ann_type()
                ann_shape = self._get_output_ann_shape()
                if ann_type is None:
                    self.logger.debug("Output annotation type: {0}".format(ann_type))
                    self.logger.debug("Output annotation shape: {0}".format(ann_shape))
                    dt = AnnotatedDataTyper(
                        v, annotated_type=ann_type, annotated_shape=ann_shape
                    )
                    t = dt.get_type()
                    self.logger.debug("Resolved annotation: {0}".format(t))
                else:
                    self.logger.debug("No annotated metadata could be retrieved")
                    t = None
                if t is None:
                    dt = PureDataTyper(v)
                    t = dt.get_type()
                output_schema[k] += [t]
        input_schema_ = {}
        for k, v in input_schema.items():
            if self.__dicts_are_identical(v):
                input_schema_[k] = v[0]
            else:
                self.logger.error("Input data types are not consistent")
        meta = self.model.autoservice._latest_meta
        self.logger.debug("Latest meta: {0}".format(meta))
        output_schema_ = {}
        for k, v in output_schema.items():
            if ann_shape == "Single":
                if self.__dicts_are_identical(v):
                    output_schema_[k] = v[0]
                else:
                    self.logger.error("Output data types are not consistent")
                    raise OutputDataTypesNotConsistentError
            elif ann_shape == "List":
                if self.__dicts_are_identical(v):
                    output_schema_[k] = v[0]
                else:
                    self.logger.warning(
                        "Output data types are not consistent. Continuing anyway, but removing output shape."
                    )
                    w = v[0]
                    if k in meta:
                        try:
                            n = len(meta[k])
                        except:
                            n = None
                    else:
                        n = None
                    w["shape"] = (n,)
                    output_schema_[k] = w
            else:
                output_schema_[k] = v[0]
        for k, v in output_schema_.items():
            self.logger.debug("{0} : {1}".format(k, v))
            meta_k = meta[k]
            self.logger.debug("Meta k: {0}".format(meta_k))
            if output_schema_[k] is None:
                output_schema_[k] = {
                    "meta": meta[k],
                    "type": None,
                }  # TODO revise if type=None is the best option
            else:
                output_schema_[k]["meta"] = meta[k]
        schema = {"input": input_schema_, "output": output_schema_}
        self.logger.debug("Schema: {0}".format(schema))
        self.logger.debug("Done with the schema!")
        return schema

    @throw_ersilia_exception()
    def _get_schema_type_for_simple_run_api_case(self):
        # read metadata
        dest_dir = self._model_path(self.model_id)
        try:
            metadata = get_metadata_from_base_dir(dest_dir)
        except FileNotFoundError:
            self.logger.debug("Metadata file not available (yet)")

        # get output type from metadata file
        output_type = metadata["Output Type"]
        if len(output_type) == 1:
            self.logger.debug("Output type is {0}".format(output_type[0]))
            output_type = output_type[0]
        elif len(output_type) == 2:
            if set(output_type) == set(["Integer", "Float"]):
                output_type = "Float"
            else:
                return None
        else:
            return None
        if output_type not in ["Float", "String"]:
            return None

        # get output shape from metadata file
        output_shape = metadata["Output Shape"]
        if output_shape not in ["Single", "List"]:
            return None

        def resolve_output_meta_in_schema(output_type, output_shape):
            if output_shape == "Single" and output_type == "Float":
                return "numeric"
            if output_shape == "Single" and output_type == "String":
                return "string"
            if output_shape == "List" and output_type == "Float":
                return "numeric_array"
            if output_shape == "List" and output_type == "String":
                return "string_array"

        output_meta_in_schema = resolve_output_meta_in_schema(output_type, output_shape)

        return output_meta_in_schema

    def _try_to_resolve_output_shape(self, meta, output_type: str):
        if output_type == "numeric_array" or output_type == "string_array":
            if type(meta) is list:
                shape = (len(meta),)
                return shape
        return None

    @throw_ersilia_exception()
    def sniff(self):
        """
        Infers the model's structure by analyzing its inputs and outputs.

        This method:
        - Calculates and saves the model size.
        - Serves the model locally.
        - Generates or retrieves sample inputs.
        - Runs the model's APIs with these inputs.
        - Collects outputs and infers input/output schemas.
        - Saves the inferred schema for future use.

        Raises
        ------
        EmptyOutputError
            If the model outputs are empty.
        OutputDataTypesNotConsistentError
            If output data types are inconsistent.
        """
        self.logger.debug("Sniffing model")
        self.logger.debug("Getting model size")
        size = self._get_size_in_mb()
        self.logger.debug("Model size is {0} MB".format(size))
        path = os.path.join(self._model_path(self.model_id), MODEL_SIZE_FILE)
        with open(path, "w") as f:
            json.dump({"size": size, "units": "MB"}, f, indent=4)
        self.logger.debug("Serving model")
        self.model.autoservice.serve()
        self.logger.debug("Iterating over APIs")
        all_schemas = {}

        is_schema_done = False
        api_names = self.model.autoservice.get_apis()

        def get_results(api_name):
            self.logger.debug("Running API: {0}".format(api_name))
            self.logger.debug(self.inputs)
            results = [
                result for result in self.model.autoservice.api(api_name, self.inputs)
            ]
            self.logger.debug(str(results)[:1000])
            for r in results:
                if not r["output"]:
                    raise EmptyOutputError(model_id=self.model_id, api_name=api_name)
            return results

        # try to get schema without making calculations, just reading metadata and output file
        # this only works when the 'run' API is the only one.
        if len(api_names) == 1 and api_names[0] == "run":
            schema_type_backup = self._get_schema_type_for_simple_run_api_case()
            self.logger.debug("This is the schema {0}".format(schema_type_backup))

        if not is_schema_done:
            for api_name in self.model.autoservice.get_apis():
                self.logger.debug("Getting schema for API {0}...".format(api_name))
                results = get_results(api_name)
                schema = self._get_schema(results)
                self.logger.debug("This is the schema {0}".format(schema))
                if api_name == "run":
                    if "outcome" in schema["output"]:
                        if schema["output"]["outcome"]["type"] is None:
                            schema["output"]["outcome"]["type"] = schema_type_backup
                        if "shape" not in schema["output"]["outcome"]:
                            shape = self._try_to_resolve_output_shape(
                                schema["output"]["outcome"]["meta"],
                                schema["output"]["outcome"]["type"],
                            )
                            if shape is not None:
                                schema["output"]["outcome"]["shape"] = shape
                all_schemas[api_name] = schema

        path = os.path.join(self._model_path(self.model_id), API_SCHEMA_FILE)
        with open(path, "w") as f:
            json.dump(all_schemas, f, indent=4)
        self.logger.debug("API schema saved at {0}".format(path))
