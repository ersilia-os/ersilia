import csv
import json
import math
import numpy as np
import os
import random
from pathlib import Path
from typing import Any

# ruff: noqa
MISSING_PACKAGES = False
try:
    from fuzzywuzzy import fuzz
    from scipy.stats import spearmanr
except ImportError:
    MISSING_PACKAGES = True
# ruff: enable
from .constants import Checks, Options, STATUS_CONFIGS, PREDEFINED_COLUMN_FILE
from .io import IOService
from ....hub.fetch.actions.template_resolver import TemplateResolver
from ....utils.exceptions_utils import test_exceptions as texc
from ....utils.hdf5 import Hdf5DataLoader
from ....utils.exceptions_utils.base_information_exceptions import _read_default_fields
from ....store.utils import echo_exceptions, ClickInterface

class CheckService:
    """
    Service for performing various checks on the model.

    Parameters
    ----------
    logger : logging.Logger
        Logger for logging messages.
    model_id : str
        Identifier of the model.
    dir : str
        Directory where the model repository is located.
    from_github : bool
        Flag indicating whether to fetch the repository from GitHub.
    from_dockerhub : bool
        Flag indicating whether to fetch the repository from DockerHub.
    ios : IOService
        Instance of IOService for handling input/output operations.

    Examples
    --------
    .. code-block:: python

        check_service = CheckService(
            logger=logger,
            model_id="model_id",
            dir="/path/to/dir",
            from_github=True,
            from_dockerhub=False,
            ios=ios,
        )
        check_service.check_files()
    """

    def __init__(
        self,
        logger: Any,
        model_id: str,
        dir: str,
        from_github: bool,
        from_dockerhub: bool,
        ios: IOService,
    ):
        self.logger = logger
        self.model_id = model_id
        self.dir = dir
        self.from_github = from_github
        self.from_dockerhub = from_dockerhub
        self.ios = ios
        self._run_check = self.ios._run_check
        self._read_metadata = self.ios._read_metadata
        self._generate_table = self.ios._generate_table
        self.get_file_requirements = self.ios.get_file_requirements
        self.console = ios.console
        self.original_smiles_list = []
        self.check_results = ios.check_results
        self.resolver = TemplateResolver(model_id=model_id, repo_path=self.dir)
        # Field defaults
        self.valid_tasks = set(_read_default_fields("Task"))
        self.valid_model_outputs = set(_read_default_fields("Output"))
        self.valid_model_inputs = set(_read_default_fields("Input"))

    def _get_output_consistency(self):
        if self.dir is not None:
            output_consistency =self.ios.get_output_consistency()
        else:
            output_consistency = "Fixed"
        return output_consistency

    def _check_file_existence(self, path):
        if not os.path.exists(os.path.join(self.dir, path)):
            raise FileNotFoundError(f"File '{path}' does not exist.")

    def check_files(self):
        """
        Check the existence of required files for the model.
        """
        requirements = self.get_file_requirements()
        for file in requirements:
            self.logger.debug(f"Checking file: {file}")
            self._run_check(self._check_file_existence, None, f"{file}", file)

    def _check_model_id(self, data):
        self.logger.debug("Checking model ID...")
        if data["Identifier"] != self.model_id:
            raise texc.WrongCardIdentifierError(self.model_id)

    def _check_model_slug(self, data):
        self.logger.debug("Checking model slug...")
        if not data["Slug"]:
            raise texc.EmptyField("slug")

    def _check_model_description(self, data):
        self.logger.debug("Checking model description...")
        if not data["Description"]:
            raise texc.EmptyField("Description")

    def _check_model_tag(self, data):
        self.logger.debug("Checking model tag...")
        if not data["Tag"]:
            raise texc.EmptyField("Tag")

    def _check_model_source_code(self, data):
        self.logger.debug("Checking model source code...")
        if not data["Source Code"]:
            raise texc.EmptyField("Source Code")

    def _check_model_source_title(self, data):
        self.logger.debug("Checking model title...")
        if not data["Title"]:
            raise texc.EmptyField("Title")

    def _check_model_status(self, data):
        self.logger.debug("Checking model status...")
        if not data["Status"]:
            raise texc.EmptyField("Status")

    def _check_model_contributor(self, data):
        key = "Contributor"
        self.logger.debug(f"Checking {key} field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)

    def _check_model_interpret(self, data):
        self.logger.debug("Checking model interpretation...")
        if not data["Interpretation"]:
            raise texc.EmptyField("Interpretation")

    def _check_model_dockerhub_url(self, data):
        key = "DockerHub"
        self.logger.debug(f"Checking {key} URL field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)

    def _check_model_s3_url(self, data):
        key = "S3"
        self.logger.debug(f"Checking {key} URL field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)

    def _check_model_arch(self, data):
        key = "Docker Architecture"
        self.logger.debug(f"Checking {key} field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)

    def _check_model_publication(self, data):
        key = "Publication"
        self.logger.debug(f"Checking {key} field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)

    def _check_model_task(self, data):
        self.logger.debug("Checking model task...")
        raw_tasks = data.get("Task")
        if isinstance(raw_tasks, str):
            tasks = [task.strip() for task in raw_tasks.split(",") if task.strip()]
        elif isinstance(raw_tasks, list):
            tasks = [
                task.strip()
                for task in raw_tasks
                if isinstance(task, str) and task.strip()
            ]
        else:
            raise texc.InvalidEntry("Task")

        if not tasks:
            raise texc.InvalidEntry("Task")

        invalid_tasks = [task for task in tasks if task not in self.valid_tasks]
        if invalid_tasks:
            raise texc.InvalidEntry("Task")

        self.logger.debug("All tasks are valid.")

    def _check_model_output(self, data):
        self.logger.debug("Checking model output...")
        raw_outputs = data.get("Output")
        if isinstance(raw_outputs, str):
            outputs = [
                output.strip() for output in raw_outputs.split(",") if output.strip()
            ]
        elif isinstance(raw_outputs, list):
            outputs = [
                output.strip()
                for output in raw_outputs
                if isinstance(output, str) and output.strip()
            ]
        else:
            raise texc.InvalidEntry("Output")

        if not outputs:
            raise texc.InvalidEntry("Output")

        invalid_outputs = [
            output for output in outputs if output not in self.valid_model_outputs
        ]
        if invalid_outputs:
            raise texc.InvalidEntry("Output")

        self.logger.debug("All outputs are valid.")

    def _check_model_input(self, data):
        self.logger.debug("Checking model input")

        model_input = data.get("Input")
        if isinstance(model_input, str):
            model_input = [
                input.strip() for input in model_input.split(",") if input.strip()
            ]
        elif isinstance(model_input, list):
            model_input = [
                input.strip()
                for input in model_input
                if isinstance(input, str) and input.strip()
            ]
        else:
            raise texc.InvalidEntry("Input")

        if not model_input:
            raise texc.InvalidEntry("Output")

        invalid_inputs = [
            input for input in model_input if input not in self.valid_model_inputs
        ]
        if invalid_inputs:
            raise texc.InvalidEntry("Input")

        self.logger.debug("All Inputs are valid.")

    def _check_model_output_type(self, data):
        self.logger.debug("Checking model output type...")
        valid_output_types = ["String", "Float", "Integer"]

        model_output_type = data.get("Output Type")
        model_output_type = (
            model_output_type[0]
            if isinstance(model_output_type, list)
            else model_output_type
        )
        if not model_output_type or model_output_type not in valid_output_types:
            raise texc.InvalidEntry("Output Type")

    # NEW FILEDS:

    def _check_model_source(self, data):
        key = "Source"
        self.logger.debug(f"Checking {key}  field..")
        if not data[key]:
            raise texc.EmptyField(key)

    def _check_model_source_type(self, data):
        key = "Source Type"
        self.logger.debug(f"Checking {key}  field..")
        if not data[key]:
            raise texc.EmptyField(key)

    def _check_model_sub_tasks(self, data):
        key = "Subtask"
        self.logger.debug(f"Checking {key}  field..")
        if not data[key]:
            raise texc.EmptyField(key)

    def _check_model_input_dim(self, data):
        key = "Input Dimension"
        self.logger.debug(f"Checking {key}  field..")
        if not data[key]:
            raise texc.EmptyField(key)

    def _check_model_target_organism(self, data):
        key = "Target Organism"
        self.logger.debug(f"Checking {key}  field..")
        if not data[key]:
            raise texc.EmptyField(key)

    def _check_model_biomedical_area(self, data):
        key = "Biomedical Area"
        self.logger.debug(f"Checking {key}  field..")
        if not data[key]:
            raise texc.EmptyField(key)

    def _check_model_output_dim(self, data):
        key = "Output Dimension"
        self.logger.debug(f"Checking {key}  field..")
        if not data[key]:
            raise texc.EmptyField(key)

    def _check_model_output_consistency(self, data):
        key = "Output Consistency"
        self.logger.debug(f"Checking {key}  field..")
        if not data[key]:
            raise texc.EmptyField(key)

    def _check_model_contribution_date(self, data):
        key = "Incorporation Date"
        self.logger.debug(f"Checking {key}  field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)

    def _check_model_image_size(self, data):
        key = "Image Size"
        self.logger.debug(f"Checking {key}  field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)

    def _check_model_env_size(self, data):
        key = "Environment Size"
        self.logger.debug(f"Checking {key}  field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)

    def _check_model_model_size(self, data):
        key = "Model Size"
        self.logger.debug(f"Checking {key}  field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)

    def _check_model_computational_performance_one(self, data):
        key = "Computational Performance 1"
        self.logger.debug(f"Checking {key}  field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)

    def _check_model_computational_performance_two(self, data):
        key = "Computational Performance 2"
        self.logger.debug(f"Checking {key}  field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)

    def _check_model_computational_performance_three(self, data):
        key = "Computational Performance 3"
        self.logger.debug(f"Checking {key}  field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)
    def _check_model_computational_performance_four(self, data):
        key = "Computational Performance 4"
        self.logger.debug(f"Checking {key}  field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)
        
    def _check_model_computational_performance_five(self, data):
        key = "Computational Performance 5"
        self.logger.debug(f"Checking {key}  field..")
        if key in data:
            if not data[key]:
                raise texc.EmptyField(key)
        else:
            raise texc.EmptyKey(key)


    def check_information(self):
        """
        Perform various checks on the model information.

        Parameters
        ----------
        output : file-like object
            The output file to write to.
        """
        self.logger.debug(f"Beginning checks for {self.model_id} model information")
        data = self._read_metadata()

        self._run_check(self._check_model_id, data, "Model ID")
        self._run_check(self._check_model_slug, data, "Model Slug")
        self._run_check(self._check_model_status, data, "Model Status")
        self._run_check(self._check_model_source_title, data, "Model Title")
        self._run_check(self._check_model_description, data, "Model Description")
        self._run_check(self._check_model_task, data, "Model Task")
        self._run_check(self._check_model_input, data, "Model Input")
        self._run_check(self._check_model_output, data, "Model Output")
        self._run_check(self._check_model_interpret, data, "Model Interpretation")
        self._run_check(self._check_model_tag, data, "Model Tag")
        self._run_check(self._check_model_publication, data, "Model Publication")
        self._run_check(self._check_model_source_code, data, "Model Source Code")
        self._run_check(self._check_model_contributor, data, "Model Contributor")
        self._run_check(self._check_model_dockerhub_url, data, "Model Dockerhub URL")
        self._run_check(self._check_model_s3_url, data, "Model S3 URL")
        self._run_check(self._check_model_arch, data, "Model Docker Architecture")
        # New added fields
        self._run_check(self._check_model_biomedical_area, data, "Model Biomodel Area")
        self._run_check(
            self._check_model_target_organism, data, "Model Target Organism"
        )
        self._run_check(self._check_model_source, data, "Model Source")
        self._run_check(self._check_model_source_type, data, "Model Source Type")
        self._run_check(self._check_model_sub_tasks, data, "Model Sub Tasks")
        self._run_check(self._check_model_input_dim, data, "Model Input Dimensions")
        self._run_check(self._check_model_output_dim, data, "Model Output Dimension")
        self._run_check(
            self._check_model_output_consistency, data, "Model Output Consistency"
        )
        self._run_check(
            self._check_model_contribution_date, data, "Model Contribution Date"
        )
        self._run_check(self._check_model_image_size, data, "Model Image Size")
        self._run_check(self._check_model_env_size, data, "Model Environment Size")
        self._run_check(self._check_model_model_size, data, "Model Directory Size")
        self._run_check(
            self._check_model_computational_performance_one,
            data,
            "Model Computational Performance for 1 input",
        )
        self._run_check(
            self._check_model_computational_performance_two,
            data,
            "Model Computational Performance for 10 input",
        )
        self._run_check(
            self._check_model_computational_performance_three,
            data,
            "Model Computational Performance for 100 input",
        )
        self._run_check(
            self._check_model_computational_performance_four,
            data,
            "Model Computational Performance for 1000 input",
        )
        self._run_check(
            self._check_model_computational_performance_five,
            data,
            "Model Computational Performance for 10000 input",
        )

    def _duplicate(self, csv_file):
        with open(csv_file, mode="r", newline="", encoding="utf-8") as file:
            reader = list(csv.DictReader(file))
            if not reader:
                return

            sr = random.choice(reader)
            duplicates = [sr.copy() for _ in range(Options.NUM_SAMPLES.value)]
        with open(csv_file, mode="a", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=reader[0].keys())
            writer.writerows(duplicates)

    def get_inputs(self, types):
        samples = IOService._get_input_from_example_file(self.dir)
        if types == "str":
            return samples[0]
        if types == "list":
            return json.dumps(samples)
        if types == "csv":
            return IOService._get_input_file_path(self.dir)

    def _is_invalid_value(self, value):
        try:
            if value is None:
                return True
            if isinstance(value, str):
                if value.strip().lower() in {"", "nan", "null", "none"}:
                    return True
            if isinstance(value, float) and (np.isnan(value) or np.isinf(value)):
                return True

            if isinstance(value, (list, tuple)):
                return any(self._is_invalid_value(item) for item in value)
            if isinstance(value, dict):
                return any(self._is_invalid_value(item) for item in value.values())
            if isinstance(value, (np.ndarray)):
                return np.any(np.isnan(value)) or np.any(np.isinf(value))
        except Exception:
            return True
        return False

    def trim_string(self, text, max_length=100):
        if len(text) <= max_length:
            return text
        return (
            text[:max_length].rsplit(" ", 1)[0] + "..."
            if " " in text[:max_length]
            else text[:max_length] + "..."
        )
    
    def _input_matchs_output_order(self, input_smiles, output_smiles):
        return input_smiles == output_smiles
    
    def _check_csv(self, file_path, input_type="list"):
        
        self.logger.debug(f"Checking CSV file: {file_path} for {input_type} input")
        error_details = []
        is_online = False
        is_fixed = False
        
        metadata = {}
        try:
            metadata = self.ios._read_metadata()
        except:
            self.logger.info("Models are running from playground!")
            metadata["Source"] = "Local"

        if "Source" in metadata:
            if metadata["Source"] == "Online":
                is_online = True
        if "Output Consistency" in metadata:
            if metadata["Output Consistency"] == "Fixed":
                is_fixed = True
        try:
            with open(file_path, "r") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                output_smiles = [row.get("input") or row.get("smiles") for row in rows]
                is_order_matchs = self._input_matchs_output_order(self.original_smiles_list, output_smiles)
                if not is_order_matchs:
                    error_details.append("Inputs and outputs order does not match!")
            if is_fixed:
                missing_values_in_first_col = self.find_missing_first_output_col(file_path)
                if not is_online:
                    if len(missing_values_in_first_col) > 0:
                        error_details.extend(missing_values_in_first_col)
                else:
                    if len(missing_values_in_first_col) == len(rows):
                        error_details.extend(missing_values_in_first_col)

            else:
                status = self._check_all_columns_not_null(file_path)
                if status[0][-1] == str(STATUS_CONFIGS.FAILED):
                    error_details.append(
                            f"All outputs are None. Invalid output content. {status[0][1]}"
                        )
                    
            if error_details:
                return (
                    f"{input_type.upper()}-CSV",
                    f"Validation failed: {', '.join(error_details)}",
                    str(STATUS_CONFIGS.FAILED),
                )
            return (
                f"{input_type.upper()}-CSV",
                "Valid Content and Input Match",
                str(STATUS_CONFIGS.PASSED),
            )
        except Exception as e:
            return (
                f"{input_type.upper()}-CSV",
                f"Validation error: {str(e)}",
                str(STATUS_CONFIGS.FAILED),
            )

    def _get_original_smiles_list(self, inp_type, inp_data):
        if inp_type == "str":
            return [inp_data]
        elif inp_type == "list":
            return json.loads(inp_data)
        elif inp_type == "csv":
            with open(inp_data, "r") as f:
                reader = csv.DictReader(f)
                return [
                    row[key]
                    for row in reader
                    for key in ("input", "smiles")
                    if key in row
                ]
        else:
            raise ValueError(f"Unsupported input type: {inp_type}")

    def check_model_output_content(self, run_example, run_model):
        status = []
        self.logger.debug("Checking model output...")
        for inp_type in Options.INPUT_TYPES_TEST.value:
            for output_file in Options.OUTPUT_FILES_TEST.value:
                inp_data = self.get_inputs(inp_type)
                self.original_smiles_list = self._get_original_smiles_list(
                    inp_type, inp_data
                )
                run_model(inputs=inp_data, output=output_file, batch=100)
                self.logger.info("validating output")
                _status = self.validate_file_content(output_file, inp_type)
                status.append(_status)
        return status

    def validate_file_content(self, file_path, input_type):
        """
        Validate output file content and check input SMILES match.
        """

        def check_json():
            self.logger.debug(
                f"Checking JSON file: {file_path} for input: {input_type}"
            )
            error_details = []
            try:
                with open(file_path, "r") as f:
                    content = json.load(f)

                def _validate_item(item, path="result"):
                    if self._is_invalid_value(item):
                        self.logger.error(f"Invalid value at {path}: {item}")
                        error_details.append(f"Invalid value at {path}: {item}")
                    elif isinstance(item, dict):
                        for key, value in item.items():
                            _validate_item(value, f"{path}.{key}")
                    elif isinstance(item, list):
                        for idx, val in enumerate(item):
                            _validate_item(val, f"{path}[{idx}]")

                _validate_item(content)

                try:
                    output_smiles = []
                    for item in content:
                        if "input" in item and isinstance(item["input"], dict):
                            output_smiles.append(item["input"].get("input", None))
                        else:
                            error_details.append(
                                "Missing 'input' structure in JSON item"
                            )
                            break

                    self.logger.info(
                        f"Matching the input SMILES in JSON file: {output_smiles} and original smiles: {self.original_smiles_list}"
                    )

                    total_outputs = len(output_smiles)
                    null_count = sum(1 for smile in output_smiles if smile is None)
                    null_percentage = (
                        (null_count / total_outputs) if total_outputs > 0 else 0
                    )

                    if self._get_output_consistency() != "Fixed":
                        if null_percentage > 0.99:
                            error_details.append(
                                "Null output percentage exceeds 25% for variable output consistency."
                            )

                        non_null_output = [s for s in output_smiles if s is not None]
                        non_null_expected = [
                            s for s in self.original_smiles_list if s is not None
                        ]
                        if non_null_output != non_null_expected:
                            error_details.append(
                                "Non-null input SMILES mismatch or order incorrect in JSON."
                            )
                    else:
                        if null_count > 0:
                            error_details.append(
                                "Null outputs found in fixed output consistency."
                            )
                        elif output_smiles != self.original_smiles_list:
                            error_details.append(
                                "Input SMILES mismatch or order incorrect in JSON."
                            )
                except Exception as e:
                    error_details.append(f"Error checking SMILES in JSON: {str(e)}")

                if error_details:
                    return (
                        f"{input_type.upper()}-JSON",
                        f"Validation failed: {', '.join(error_details)}",
                        str(STATUS_CONFIGS.FAILED),
                    )
                return (
                    f"{input_type.upper()}-JSON",
                    "Valid Content and Input Match",
                    str(STATUS_CONFIGS.PASSED),
                )
            except Exception as e:
                return (
                    f"{input_type.upper()}-JSON",
                    f"Validation error: {str(e)}",
                    str(STATUS_CONFIGS.FAILED),
                )

        def check_h5():
            
            self.logger.debug(f"Checking HDF5 file: {file_path}")
            error_details = []
            metadata = {}
            try:
                metadata = self.ios._read_metadata()
            except:
                self.logger.info("Models are running from playground!")
                metadata["Source"] = "Local"

            is_online = False
            if "Source" in metadata:
                if metadata["Source"] == "Online":
                    is_online = True
            is_fixed = False
            if "Output Consistency" in metadata:
                if metadata["Output Consistency"] == "Fixed":
                    is_fixed = True
            if is_fixed and is_online:
                return (
                    f"{input_type.upper()}-HDF5",
                    f"Test not applicable for online models with fixed output consistency",
                    str(STATUS_CONFIGS.WARNING),
                )
            
            try:
                loader = Hdf5DataLoader()
                loader.load(file_path)

                content = next(
                    (
                        x
                        for x in [
                            loader.values,
                            loader.keys,
                            loader.inputs,
                            loader.features,
                        ]
                        if x is not None
                    ),
                    None,
                )
                
                output_smiles = (
                    [s for s in loader.inputs] if loader.inputs is not None else []
                )
                is_order_matchs = self._input_matchs_output_order(self.original_smiles_list, output_smiles)
                
                if not is_order_matchs:
                    error_details.append("Order Mis match happens")

                if content is None or (hasattr(content, "size") and content.size == 0):
                    error_details.append("Empty content")

                content_array = np.array(content)
                    

                if np.issubdtype(content_array.dtype, np.floating):
                    invalid_mask = np.isnan(content_array)
                elif np.issubdtype(content_array.dtype, np.integer):
                    invalid_mask = np.full(content_array.shape, False)
                elif content_array.dtype.kind in ['S', 'U']:
                    invalid_mask = (content_array == '')
                elif content_array.dtype.kind == 'O':
                    def is_empty(x):
                        if isinstance(x, bytes):
                            return x == b''
                        elif isinstance(x, str):
                            return x == ''
                        else:
                            return False

                    vector_is_empty = np.vectorize(is_empty)
                    invalid_mask = vector_is_empty(content_array)
                else:
                    invalid_mask = np.full(content_array.shape, False)

                if invalid_mask.any():
                    invalid_indices = np.argwhere(invalid_mask)
                    error = "H5 content invalid value at index: " + ", ".join(str(index) for index in invalid_indices)
                    error = self.trim_string(error)
                    self.logger.error(error)
                    error_details.append(error)

                if error_details:
                    return (
                        f"{input_type.upper()}-HDF5",
                        error_details,
                        str(STATUS_CONFIGS.FAILED),
                    )

                return (
                    f"{input_type.upper()}-HDF5",
                    "Valid content and Input Match"
                    if not error_details
                    else f"Errors: {', '.join(error_details)}",
                    str(
                        (STATUS_CONFIGS.PASSED)
                        if not error_details
                        else STATUS_CONFIGS.FAILED
                    ),
                )

            except Exception as e:
                return (
                    f"{input_type.upper()}-HDF5",
                    f"Validation error: {str(e)}",
                    str(STATUS_CONFIGS.FAILED),
                )

        if not Path(file_path).exists():
            raise FileNotFoundError(f"File {file_path} not found.")

        file_ext = Path(file_path).suffix.lower()
        if file_ext == ".json":
            return check_json()
        elif file_ext == ".h5":
            return check_h5()
        elif file_ext == ".csv":
            return self._check_csv(file_path, input_type)
        else:
            raise ValueError(f"Unsupported file type: {file_ext}")

    def _read_column_header(self, reader):
        return [row[0] for row in reader if row][1:]
    
    def _find_csv_mismatches(self, csv_out_one, csv_out_two):        
        is_online = False
        metadata = {}
        try:
            metadata = self.ios._read_metadata()
        except:
            self.logger.info("Models are running from playground!")
            metadata["Source"] = "Local"
        if "Source" in metadata:
            if metadata["Source"] == "Online":
                is_online = True

        with open(csv_out_one) as f:
            reader = csv.reader(f)
            header = next(reader)
            idxs = [i for i, col in enumerate(header) if col != "input" and col != "key"]
            rows1 = []
            for r in reader:
                if len(r) != len(header):
                    raise Exception("There was a row with less columns than expected")
                rows1.append([r[i] for i in idxs])
        with open(csv_out_two) as f:
            reader = csv.reader(f)
            header = next(reader)
            idxs = [i for i, col in enumerate(header) if col != "input" and col != "key"]
            rows2 = []
            for r in reader:
                if len(r) != len(header):
                    return [(Checks.COLUMN_MISMATCH, "There was a row with less columns than expected", str(STATUS_CONFIGS.FAILED))]
                rows2.append([r[i] for i in idxs])
        mismatches = []
        max_rows = max(len(rows1), len(rows2))
        for i in range(max_rows):
            row1 = rows1[i] if i < len(rows1) else []
            row2 = rows2[i] if i < len(rows2) else []
            max_cols = max(len(row1), len(row2))
            is_empty_row1 = True
            is_empty_row2 = True
            for x in row1:
                if x:
                    is_empty_row1 = False
                    break
            for x in row2:
                if x:
                    is_empty_row2 = False
                    break
            if is_empty_row1 or is_empty_row2:
                if is_online:
                    continue
            for j in range(max_cols):
                v1 = row1[j] if j < len(row1) else None
                v2 = row2[j] if j < len(row2) else None

                try:
                    if v1 == "" or v2 == "":
                        raise ValueError

                    f1 = float(v1)
                    f2 = float(v2)

                    if math.isnan(f1) and math.isnan(f2):
                        continue

                    if abs(f1-f2)>0.1:
                        mismatches.append((i, j, v1, v2))

                except (ValueError, TypeError):
                    if v1 != v2:
                        mismatches.append((i, j, v1, v2))

        if mismatches:
            if len(mismatches) > 3:
                return [(Checks.COLUMN_MISMATCH, f"Column mismatches found (and more): {mismatches[:3]}", str(STATUS_CONFIGS.FAILED))]
            else:
                return [(Checks.COLUMN_MISMATCH, f"Column mismatch found: {mismatches}", str(STATUS_CONFIGS.FAILED))]
        return [(Checks.COLUMN_MISMATCH, "No column mismatches", str(STATUS_CONFIGS.PASSED))]
    
    def find_missing_first_output_col(self, path):
        missing = []
        with open(path, newline='') as f:
            reader = csv.reader(f)
            for i, row in enumerate(reader):
                val = row[2] 
                if self._is_invalid_value(val):
                    missing.append((i, 3))
        return missing
    
    def _check_all_columns_not_null(self, path):
        col_index = 2
        missing_rows = []
        total_rows = 0
        with open(path, newline='', encoding='utf-8') as f:
            reader = csv.reader(f)
            missing_cols = []
            for i, row in enumerate(reader):
                total_rows += 1
                vals = row[col_index:]
                num_missing = 0
                for val in vals:                
                    if self._is_invalid_value(val):
                        num_missing += 1
                missing_cols += [num_missing]
            total_cols = len(row) - 2

        fully_empty_rows = [i for i, val in enumerate(missing_cols) if val == total_cols]

        partially_empty_rows = [i for i, val in enumerate(missing_cols) if val > 0 and val < total_cols]

        if len(fully_empty_rows) == total_rows:
            return [(Checks.EMPTY_COLUMNS, "All columns values are nulls", str(STATUS_CONFIGS.FAILED))]
        
        if len(fully_empty_rows) > 0:
            missing_rows += fully_empty_rows
            return [(Checks.EMPTY_COLUMNS, f"All columns values are nulls at these row indices: {missing_rows}", str(STATUS_CONFIGS.WARNING))]

        if len(partially_empty_rows) > 0:
            missing_rows += partially_empty_rows
            return [(Checks.EMPTY_COLUMNS, f"Some columns values are nulls at these row indices: {missing_rows}", str(STATUS_CONFIGS.WARNING))]
                
        return [(Checks.EMPTY_COLUMNS, "All columns values are not nulls", str(STATUS_CONFIGS.PASSED))]
            
    def compare_csv_columns(self, column_csv, csv_file):
        try:
            with open(column_csv, "r", newline="") as f1, open(csv_file, "r", newline="") as f2: # ruff: noqa: E501
                reader1 = csv.reader(f1)
                reader2 = csv.reader(f2)

                first_column_values = self._read_column_header(reader1)
                header2 = next(reader2, None)
                header2 = header2 if "framework/examples" in csv_file else header2[2:]

                if not first_column_values or header2 is None:
                    return [
                        (
                            Checks.COLUMN_NAME_VALIDITY.value,
                            Checks.COLUMN_CHECK_FAILURE.value,
                            str(STATUS_CONFIGS.FAILED),
                        )
                    ]

                if first_column_values == header2:
                    return [
                        (
                            Checks.COLUMN_NAME_VALIDITY.value,
                            Checks.COLUMN_CHECK_SUCCESS.value,
                            str(STATUS_CONFIGS.PASSED),
                        )
                    ]
                else:
                    return [
                        (
                            Checks.COLUMN_NAME_VALIDITY.value,
                            Checks.COLUMN_CHECK_FAILURE.value,
                            str(STATUS_CONFIGS.FAILED),
                        )
                    ]

        except Exception as e:
            return [
                (
                    Checks.COLUMN_NAME_VALIDITY.value,
                    f"An error occurred: {e}",
                    str(STATUS_CONFIGS.FAILED),
                )
            ]

    def check_simple_model_output(self, run_model):
        input_path = IOService._get_input_file_path(self.dir)
        output_path = IOService._get_output_file_path(self.dir)
        run_model(inputs=input_path, output=Options.OUTPUT_CSV.value, batch=100)
        output_consistency = self._get_output_consistency()
        if output_consistency == "Fixed":
            res_one = self._find_csv_mismatches(output_path, Options.OUTPUT_CSV.value)
        else:
            res_one = self._check_all_columns_not_null(Options.OUTPUT_CSV.value)
        res_two = self.compare_csv_columns(
            os.path.join(self.dir, PREDEFINED_COLUMN_FILE), Options.OUTPUT_CSV.value
        )
        _completed_status = []
        if res_one[0][-1] == str(STATUS_CONFIGS.FAILED):
            self.logger.error("Model output has content problem")
            _completed_status.append(
                (
                    Checks.SIMPLE_MODEL_RUN.value,
                    res_one[0][1],
                    str(STATUS_CONFIGS.FAILED),
                )
            )
            return _completed_status
        _completed_status.append(
            (
                Checks.SIMPLE_MODEL_RUN.value,
                res_one[0][1],
                str(STATUS_CONFIGS.PASSED),
            )
        )
        first_item = res_two[0]
        _res_two = (Checks.SIMPLE_MODEL_RUN_COLUMNS.value,) + first_item[1:]
        res_two[0] = _res_two
        _completed_status.extend(res_two)
        return _completed_status

    
    def check_consistent_output(self, run_example, run_model):
        """
        Check if the model produces consistent output.

        Parameters
        ----------
        run_example : callable
            Function to generate example input.
        run_model : callable
            Function to run the model.
        """
        self.logger.debug("Confirming model produces consistent output...")

        def compute_rmse(y_true, y_pred):
            return sum((yt - yp) ** 2 for yt, yp in zip(y_true, y_pred)) ** 0.5 / len(
                y_true
            )

        def _compare_output_strings(output1, output2):
            return fuzz.ratio(output1, output2)

        def validate_output(output1, output2):
            if isinstance(output1, (float, int)):
                rmse = compute_rmse([output1], [output2])
                if rmse > 0.1:
                    echo_exceptions(f"Model output is inconsistent. The RMSE between two consequetive outputs found to be > 10%:[{rmse*100}]", ClickInterface())
                    raise ValueError

                rho, _ = spearmanr([output1], [output2])
                if rho < 0.5:
                    echo_exceptions(f"Model output is inconsistent. The Spearman correlation between two consecutive outputs found to be < 50%", ClickInterface())
                    raise ValueError

            elif isinstance(output1, list):
                rmse = compute_rmse(output1, output2)
                if rmse > 0.1:
                    echo_exceptions(f"Model output is inconsistent. The RMSE between two consecutive outputs found to be > 10%:[{rmse*100}]", ClickInterface())
                    raise ValueError
                
                rho, _ = spearmanr(output1, output2)

                if rho < 0.5:
                    echo_exceptions(f"Model output is inconsistent. The Spearman correlation between two consecutive outputs found to be < 50%", ClickInterface())
                    raise ValueError
                
            elif isinstance(output1, str):
                if _compare_output_strings(output1, output2) <= 95:
                    echo_exceptions(f"Model output is inconsistent. The Fuzz ratio correlation between two consecutive outputs found to be <= 95%", ClickInterface())
                    raise ValueError
                
        def read_csv(file_path):
            absolute_path = os.path.abspath(file_path)
            if not os.path.exists(absolute_path):
                echo_exceptions(f"File not found: {absolute_path}", ClickInterface())

            with open(absolute_path, mode="r") as csv_file:
                reader = csv.DictReader(csv_file)
                self.logger.info(f"Reading csv at consistency outout section: {reader}")
                return [row for row in reader]

        output1_path = os.path.abspath(Options.OUTPUT1_CSV.value)
        output2_path = os.path.abspath(Options.OUTPUT2_CSV.value)

        self.logger.debug("Confirming model produces consistent output...")

        input = self.get_inputs(types="csv")

        run_model(inputs=input, output=output1_path, batch=100)
        run_model(inputs=input, output=output2_path, batch=100)
        self.original_smiles_list = self._get_original_smiles_list("csv", input)
        check_status_one = self._check_csv(output1_path, input_type="csv")
        self.original_smiles_list = self._get_original_smiles_list("csv", input)
        check_status_two = self._check_csv(output2_path, input_type="csv")
        _completed_status = []
        if check_status_one[-1] == str(STATUS_CONFIGS.FAILED) or check_status_two[
            -1
        ] == str(STATUS_CONFIGS.FAILED):
            self.logger.error("Model output has content problem")
            _completed_status.append(
                (
                    Checks.MODEL_CONSISTENCY.value,
                    check_status_one[1],
                    str(STATUS_CONFIGS.FAILED),
                )
            )
            return _completed_status
        else:
            data1 = read_csv(output1_path)
            data2 = read_csv(output1_path)

            try:
                for res1, res2 in zip(data1, data2):
                    for key in res1:
                        if key in res2:
                            validate_output(res1[key], res2[key])
                        else:
                            raise KeyError(f"Key '{key}' not found in second result.")
                _completed_status.append(
                    (
                        Checks.MODEL_CONSISTENCY.value,
                        Checks.CONSISTENCY.value,
                        str(STATUS_CONFIGS.PASSED),
                    )
                )
            except ValueError:
                return _completed_status.append(
                    (
                        Checks.MODEL_CONSISTENCY.value,
                        Checks.INCONSISTENCY.value,
                        str(STATUS_CONFIGS.FAILED),
                    )
                )
        return _completed_status
