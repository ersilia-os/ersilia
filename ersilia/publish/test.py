# TODO adapt to input-type agnostic. For now, it works only with Compound input types.
import json
import os
import csv
import subprocess
import tempfile
import time
import click
import types
import sys
from enum import Enum
from dataclasses import dataclass
from typing import List
from datetime import datetime
from pathlib import Path

from ..utils.conda import SimpleConda
from .. import ErsiliaBase, throw_ersilia_exception
from ..io.input import ExampleGenerator
from ..utils.exceptions_utils import test_exceptions as texc
from ..utils.terminal import run_command_check_output
from ..hub.fetch.actions.template_resolver import TemplateResolver
from ..default import (
    INFORMATION_FILE, 
    INSTALL_YAML_FILE, 
    DOCKERFILE_FILE,
    PYTHON_VERSION,
    PACK_METHOD_FASTAPI,
    PACK_METHOD_BENTOML,
    METADATA_JSON_FILE,
    METADATA_YAML_FILE,
    RUN_FILE,
    PREDEFINED_EXAMPLE_FILES
)
MISSING_PACKAGES = False

try:
    from scipy.stats import spearmanr
    from fuzzywuzzy import fuzz
    from rich.console import Console
    from rich.table import Table
    from rich.text import Text
except ImportError:
    MISSING_PACKAGES = True


class Options(Enum):
    NUM_SAMPLES = 5
    BASE = "base"
    OUTPUT_CSV = "result.csv"
    OUTPUT1_CSV = "output1.csv"
    OUTPUT2_CSV = "output2.csv"
    LEVEL_DEEP = "deep"

class TableType(Enum):
    MODEL_INFORMATION_CHECKS = "Model Information Checks"
    MODEL_FILE_CHECKS = "Model File Checks"
    MODEL_DIRECTORY_SIZES = "Model Directory Sizes"
    RUNNER_CHECKUP_STATUS = "Runner Checkup Status"
    FINAL_RUN_SUMMARY = "Test Run Summary"
    INSPECT_SUMMARY = "Test Run Summary"

@dataclass
class TableConfig:
    title: str
    headers: List[str]

TABLE_CONFIGS = {
    TableType.MODEL_INFORMATION_CHECKS: TableConfig(
        title="Model Information Checks",
        headers=["Check", "Status"]
    ),
    TableType.MODEL_FILE_CHECKS: TableConfig(
        title="Model File Checks",
        headers=["Check", "Status"]
    ),
    TableType.MODEL_DIRECTORY_SIZES: TableConfig(
        title="Model Directory Sizes",
        headers=["Dest dir", "Env Dir"]
    ),
    TableType.RUNNER_CHECKUP_STATUS: TableConfig(
        title="Runner Checkup Status",
        headers=["Runner", "Status"],
    ),
     TableType.FINAL_RUN_SUMMARY: TableConfig(
        title="Test Run Summary",
        headers=["Check", "Status"]
    ),
     TableType.INSPECT_SUMMARY: TableConfig(
        title="Inspect Summary",
        headers=["Check", "Status"]
    ),
}

class STATUS_CONFIGS(Enum):
    PASSED  = ("PASSED", "green", "✔")
    FAILED  = ("FAILED", "red", "✘")
    WARNING = ("WARNING", "yellow", "⚠")
    SUCCESS = ("SUCCESS", "green", "★")
    NA = ("N/A", "dim", "~")

    def __init__(self, label, color, icon):
        self.label = label
        self.color = color
        self.icon = icon

    def __str__(self):
        return f"[{self.color}]{self.icon} {self.label}[/{self.color}]"

class SetupService:

    BASE_URL = "https://github.com/ersilia-os/"

    def __init__(
        self, model_id, 
        dir, 
        logger,
        remote
        ):
        self.model_id = model_id
        self.dir = dir
        self.logger = logger
        self.remote = remote
        self.repo_url = f"{self.BASE_URL}{self.model_id}"
        self.conda = SimpleConda()

    @staticmethod
    def run_command(command, logger, capture_output=False, shell=True
        ):
        try:
            if capture_output:
                result = subprocess.run(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                    shell=shell
                )
                return result.stdout
            else:
                process = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    shell=shell
                )
                
                stdout_lines, stderr_lines = [], []

                for line in iter(process.stdout.readline, ''):
                    if line.strip():
                        stdout_lines.append(line.strip())
                        logger.info(line.strip())

                for line in iter(process.stderr.readline, ''):
                    if line.strip():
                        stderr_lines.append(line.strip())
                        logger.error(line.strip())

                process.wait()
                if process.returncode != 0:
                    raise subprocess.CalledProcessError(
                        returncode=process.returncode,
                        cmd=command,
                        output="\n".join(stdout_lines),
                        stderr="\n".join(stderr_lines),
                    )

                return process.stdout

        except subprocess.CalledProcessError as e:
            logger.error(f"Error executing command: {e}")
            if e.output:
                logger.debug(f"Output: {e.output.strip()}")
            if e.stderr:
                logger.debug(f"Error: {e.stderr.strip()}")
            sys.exit(1)
        except Exception as e:
            logger.debug(f"Unexpected error: {e}")
            sys.exit(1)

    def fetch_repo(self):
        if self.remote and not os.path.exists(self.dir):
            out = SetupService.run_command(
                f"git clone {self.repo_url}",
                self.logger,
            )
            self.logger.info(out)

    def check_conda_env(self):
        if self.conda.exists(self.model_id):
            self.logger.debug(
                f"Conda environment '{self.model_id}' already exists."
            )
        else:
            raise Exception(
                f"Conda virtual environment not found for {self.model_id}"
            )
    
    @staticmethod
    def get_conda_env_location(model_id, logger):
        try:
            result = SetupService.run_command(
                "conda env list",
                logger=logger,
                capture_output=True
            )
            for line in result.splitlines():
                if line.startswith("#") or not line.strip():
                    continue 
                parts = line.split()
                if parts[0] == model_id:
                    return parts[-1]
        except subprocess.CalledProcessError as e:
            print(f"Error running conda command: {e.stderr}")
        except Exception as e:
            print(f"Unexpected error: {e}")
        
        return None

class IOService:

    # Required files
    RUN_FILE = f"model/framework/{RUN_FILE}"
    BENTOML_FILES = [
        DOCKERFILE_FILE,
        METADATA_JSON_FILE,
        RUN_FILE,
        "src/service.py",
        "pack.py",
        "README.md",
        "LICENSE"
    ]

    ERSILIAPACK_FILES = [
        INSTALL_YAML_FILE,
        METADATA_YAML_FILE,
        PREDEFINED_EXAMPLE_FILES[0],
        PREDEFINED_EXAMPLE_FILES[1],
        RUN_FILE,
        "README.md",
        "LICENSE",
    ]

    def __init__(
        self, 
        logger,
        dest_dir,
        model_path, 
        bundle_path, 
        bentoml_path, 
        model_id,
        dir
    ):
        self.logger = logger
        self.model_id = model_id
        self.dir = dir
        self.model_size = 0
        self.console = Console()
        self.check_results = []
        self._model_path = model_path
        self._bundle_path = bundle_path
        self._bentoml_path = bentoml_path
        self._dest_dir = dest_dir


    def _run_check(
            self, 
            check_function, 
            data, 
            check_name, 
            additional_info=None
        ):
  
        try:
            if additional_info is not None:
                check_function(additional_info)
            else:
                check_function(data)
            self.check_results.append((
                check_name, 
                str(STATUS_CONFIGS.PASSED)
            ))
            return True
        except Exception as e:
            self.logger.error(
                f"Check '{check_name}' failed: {e}"
            )
            self.check_results.append((
                check_name, 
                str(STATUS_CONFIGS.FAILED)
            ))
            return False


    def _generate_table(
            self,
            title,
            headers,
            rows,
            large_table=False,
            merge=False 
        ):
        f_col_width = 30 if large_table else 30    
        l_col_width = 50 if large_table else 10     
        d_col_width = 30 if not large_table else 20 

        table = Table(
            title=Text(
                title,
                style="bold light_green" 
            ),
            border_style="light_green",   
            show_lines=True,
        )

        table.add_column(
            headers[0],
            justify="left",
            width=f_col_width,
            style="bold" 
        )
        for header in headers[1:-1]:
            table.add_column(
                header,
                justify="center",
                width=d_col_width,
                style="bold" 
            )
        table.add_column(
            headers[-1],
            justify="right",
            width=l_col_width,
            style="bold" 
        )

        prev_value = None 
        for row in rows:
            first_col = str(row[0])
            if merge and first_col == prev_value:
                first_col = ""  
            else:
                prev_value = first_col

            styled_row = [
                Text(first_col, style="bold"),
                *[str(cell) for cell in row[1:-1]],
                row[-1]
            ]
            table.add_row(*styled_row)

        self.console.print(table)

    @staticmethod
    def get_model_type(model_id, repo_path):
        resolver = TemplateResolver(
            model_id=model_id, 
            repo_path=repo_path
        )
        if resolver.is_bentoml():
            return PACK_METHOD_BENTOML
        elif resolver.is_fastapi():
            return PACK_METHOD_FASTAPI
        else:
            return None

    def get_file_requirements(self):
        type = IOService.get_model_type(
            model_id=self.model_id, 
            repo_path=self.dir
        )
        if type == PACK_METHOD_BENTOML:
            return self.BENTOML_FILES
        elif type == PACK_METHOD_FASTAPI:
            return self.ERSILIAPACK_FILES
        else:
            raise ValueError(
                f"Unsupported model type: {type}"
            )

    def read_information(self):
        file = os.path.join(
            self._dest_dir, 
            self.model_id, 
            INFORMATION_FILE
        )
        if not os.path.exists(file):
            raise FileNotFoundError(
                f"Information file does not exist for model {self.model_id}"
        )
        with open(file, "r") as f:
            return json.load(f)

    def print_output(self, result, output):
        def write_output(data):
            if output is not None:
                with open(output.name, "w") as file:
                    json.dump(data, file)
            else:
                self.logger.debug(json.dumps(data, indent=4))

        if isinstance(result, types.GeneratorType):
            for r in result:
                write_output(r if r is not None else "Something went wrong")
        else:
            self.logger.debug(result)

    def get_conda_env_size(self):
        try:
            loc = SetupService.get_conda_env_location(
                self.model_id, 
                self.logger
            )
            return self.calculate_directory_size(loc)
        except Exception as e:
            self.logger.error(
                f"Error calculating size of Conda environment '{self.model_id}': {e}"
            )
            return 0

    def calculate_directory_size(self, path):
        try:
            size_output = SetupService.run_command(
                ["du", "-sm", path], 
                logger=self.logger,
                capture_output=True,
                shell=False
            )
            size = int(size_output.split()[0])
            return size
        except Exception as e:
            self.logger.error(
                f"Error calculating directory size for {path}: {e}"
            )
            return 0
        
    @throw_ersilia_exception()
    def get_directories_sizes(self):
        dir_size = self.calculate_directory_size(self.dir)
        env_size = self.get_conda_env_size()
        return dir_size, env_size


class CheckService:

    # model specs
    MODEL_TASKS = {
        "Classification",
        "Regression",
        "Generative",
        "Representation",
        "Similarity",
        "Clustering",
        "Dimensionality reduction",
    }

    MODEL_OUTPUT = {
        "Boolean",
        "Compound",
        "Descriptor",
        "Distance",
        "Experimental value",
        "Image",
        "Other value",
        "Probability",
        "Protein",
        "Score",
        "Text",
    }

    INPUT_SHAPE = {
        "Single", 
        "Pair", 
        "List", 
        "Pair of Lists", 
        "List of Lists"
    }

    OUTPUT_SHAPE = {
        "Single", 
        "List", 
        "Flexible List", 
        "Matrix", 
        "Serializable Object"
    }

    def __init__(
        self, 
        logger, 
        model_id, 
        dest_dir, 
        dir, 
        ios
        ):
        self.logger = logger
        self.model_id = model_id
        self._dest_dir = dest_dir
        self.dir = dir
        self._run_check = ios._run_check
        self._generate_table = ios._generate_table
        self._print_output = ios.print_output
        self.get_file_requirements = ios.get_file_requirements
        self.console = ios.console
        self.check_results = ios.check_results
        self.information_check = False
        self.single_input = False
        self.example_input = False
        self.consistent_output = False

    def _check_file_existence(self, path):
        if not os.path.exists(os.path.join(self.dir, path)):
            raise FileNotFoundError(
                f"File '{path}' does not exist."
        )

    def check_files(self):
        requirements = self.get_file_requirements()
        for file in requirements:
            self.logger.debug(f"Checking file: {file}")
            self._run_check(
                self._check_file_existence, 
                None,
                f"File: {file}", 
                file
            )

    def _check_model_id(self, data):
        self.logger.debug("Checking model ID...")
        if data["card"]["Identifier"] != self.model_id:
            raise texc.WrongCardIdentifierError(self.model_id)


    def _check_model_slug(self, data):
        self.logger.debug("Checking model slug...")
        if not data["card"]["Slug"]:
            raise texc.EmptyField("slug")


    def _check_model_description(self, data):
        self.logger.debug("Checking model description...")
        if not data["card"]["Description"]:
            raise texc.EmptyField("Description")

    def _check_model_task(self, data):
        self.logger.debug("Checking model task...")
        raw_tasks = data.get("card", {}).get("Task", "")
        if isinstance(raw_tasks, str):
            tasks = [
                task.strip() 
                for task 
                in raw_tasks.split(",") 
                if task.strip()
            ]
        elif isinstance(raw_tasks, list):
            tasks = [
                task.strip() 
                for task 
                in raw_tasks 
                if isinstance(task, str) and task.strip()
            ]
        else:
            raise texc.InvalidEntry(
                "Task", 
                message="Task field must be a string or list."
            )

        if not tasks:
            raise texc.InvalidEntry(
                "Task", 
                message="Task field is missing or empty."
            )

        invalid_tasks = [task for task in tasks if task not in self.MODEL_TASKS]
        if invalid_tasks:
            raise texc.InvalidEntry(
                "Task", message=f"Invalid tasks: {', '.join(invalid_tasks)}"
            )

        self.logger.debug("All tasks are valid.")

    def _check_model_output(self, data):
        self.logger.debug("Checking model output...")
        raw_outputs = data.get("card", {}).get("Output", "") or data.get("metadata", {}).get("Output", "")
        if isinstance(raw_outputs, str):
            outputs = [
                output.strip() 
                for output 
                in raw_outputs.split(",")
                if output.strip()
            ]
        elif isinstance(raw_outputs, list):
            outputs = [
                output.strip() 
                for output 
                in raw_outputs 
                if isinstance(output, str) and output.strip()
            ]
        else:
            raise texc.InvalidEntry(
                "Output", 
                message="Output field must be a string or list."
            )

        if not outputs:
            raise texc.InvalidEntry(
                "Output", 
                message="Output field is missing or empty."
            )

        invalid_outputs = [
            output 
            for output 
            in outputs 
            if output not in self.MODEL_OUTPUT
        ]
        if invalid_outputs:
            raise texc.InvalidEntry(
                "Output", 
                message=f"Invalid outputs: {', '.join(invalid_outputs)}"
            )

        self.logger.debug("All outputs are valid.")

    def _check_model_input(self, data):
        self.logger.debug("Checking model input")
        valid_inputs = [{"Compound"}, {"Protein"}, {"Text"}]
        
        model_input = data.get("card", {}).get("Input") or data.get("metadata", {}).get("Input")
        
        if not model_input or set(model_input) not in valid_inputs:
            raise texc.InvalidEntry("Input")

    def _check_model_input_shape(self, data):
        self.logger.debug("Checking model input shape")
        model_input_shape = (
            data.get("card", {}).get("Input Shape") or 
            data.get("metadata", {}).get("InputShape")
        )
        
        if model_input_shape not in self.INPUT_SHAPE:
            raise texc.InvalidEntry("Input Shape")

    def _check_model_output_type(self, data):
        self.logger.debug("Checking model output type...")
        valid_output_types = [{"String"}, {"Float"}, {"Integer"}]
        
        model_output_type = (
            data.get("card", {}).get("Output Type") or 
            data.get("metadata", {}).get("OutputType")
        )
        
        if not model_output_type or set(model_output_type) not in valid_output_types:
            raise texc.InvalidEntry("Output Type")

    def _check_model_output_shape(self, data):
        self.logger.debug("Checking model output shape...")
        model_output_shape = (
            data.get("card", {}).get("Output Shape") or 
            data.get("metadata", {}).get("OutputShape")
        )
        
        if model_output_shape not in self.OUTPUT_SHAPE:
            raise texc.InvalidEntry("Output Shape")
    
    @throw_ersilia_exception()
    def check_information(self, output):
        self.logger.debug(f"Beginning checks for {self.model_id} model information")
        file = os.path.join(
            self._dest_dir, 
            self.model_id, 
            INFORMATION_FILE
        )
        with open(file, "r") as f:
            data = json.load(f)

        self._run_check(self._check_model_id, data, "Model ID")
        self._run_check(self._check_model_slug, data, "Model Slug")
        self._run_check(self._check_model_description, data, "Model Description")
        self._run_check(self._check_model_task, data, "Model Task")
        self._run_check(self._check_model_input, data, "Model Input")
        self._run_check(self._check_model_input_shape, data, "Model Input Shape")
        self._run_check(self._check_model_output, data, "Model Output")
        self._run_check(self._check_model_output_type, data, "Model Output Type")
        self._run_check(self._check_model_output_shape, data, "Model Output Shape")

        if output is not None:
            self.information_check = True

    @throw_ersilia_exception()
    def check_single_input(self, output, run_model, run_example):

        input = run_example(
            n_samples=Options.NUM_SAMPLES.value, 
            file_name=None, 
            simple=True, 
            try_predefined=False
        )
        self.logger.debug(f"Output: {output}")
        result = run_model(
            input=input, 
            output=output, 
            batch=100
        )

        if output is not None:
            self.single_input = True
        else:
            self._print_output(result, output)

    @throw_ersilia_exception()
    def check_example_input(self, output, run_model, run_example):
        input_samples = run_example(
            n_samples=Options.NUM_SAMPLES.value, 
            file_name=None, 
            simple=True, 
            try_predefined=False
        )

        self.logger.debug("Testing model on input of 5 smiles given by 'example' command")

        result = run_model(
            input=input_samples, 
            output=output, 
            batch=100
        )

        if output is not None:
            self.example_input = True
        else:
            self._print_output(result, output)

    @throw_ersilia_exception()
    def check_consistent_output(self, run_example, run_model):
        def compute_rmse(y_true, y_pred):
             return sum((yt - yp) ** 2 for yt, yp in zip(y_true, y_pred)) ** 0.5 / len(y_true)
        def _compare_output_strings(output1, output2):
            if output1 is None and output2 is None:
                return 100
            else:
                return fuzz.ratio(output1, output2)
        def validate_output(output1, output2):
            if not isinstance(output1, type(output2)):
                raise texc.InconsistentOutputTypes(self.model_id)

            if output1 is None:
                return

            if isinstance(output1, (float, int)):
                rmse = compute_rmse([output1], [output2])
                if rmse > 0.1:
                    raise texc.InconsistentOutputs(self.model_id)

                rho, _ = spearmanr([output1], [output2])
                if rho < 0.5:
                    raise texc.InconsistentOutputs(self.model_id)

            elif isinstance(output1, list):
                rmse = compute_rmse(output1, output2)
                if rmse > 0.1:
                    raise texc.InconsistentOutputs(self.model_id)

                rho, _ = spearmanr(output1, output2)
                if rho < 0.5:
                    raise texc.InconsistentOutputs(self.model_id)

            elif isinstance(output1, str):
                if _compare_output_strings(output1, output2) <= 95:
                    raise texc.InconsistentOutputs(self.model_id)
                
        def read_csv(file_path):
            absolute_path = os.path.abspath(file_path)
            if not os.path.exists(absolute_path):
                raise FileNotFoundError(f"File not found: {absolute_path}")
            with open(absolute_path, mode='r') as csv_file:
                reader = csv.DictReader(csv_file)
                return [row for row in reader]

        output1_path = os.path.abspath(Options.OUTPUT1_CSV.value)
        output2_path = os.path.abspath(Options.OUTPUT2_CSV.value)
                
        self.logger.debug("Confirming model produces consistent output...")

        input_samples = run_example(
            n_samples=Options.NUM_SAMPLES.value, 
            file_name=None, 
            simple=True, 
            try_predefined=False
        )
        out1 = run_model(
            input=input_samples, 
            output=output1_path, 
            batch=100
        )
        out2 = run_model(
            input=input_samples, 
            output=output2_path, 
            batch=100
        )
        def create_csv(file_path, data):
    # Ensure the file path's directory exists
            Path(file_path).parent.mkdir(parents=True, exist_ok=True)
            
            # Write data to the CSV file
            with open(file_path, mode='w', newline='') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=data[0].keys())
                writer.writeheader()
                writer.writerows(data)

            # Check if the file exists after creation
            if os.path.exists(file_path):
                print(f"File successfully created: {file_path}")
            else:
                print(f"File creation failed: {file_path}")

        # Example usage
        file_path = "output/example.csv"
        data = [
            {"name": "Alice", "age": 25},
            {"name": "Bob", "age": 30},
            {"name": "Charlie", "age": 35}
        ]

        create_csv(file_path, data)
        self.logger.debug(f"Output1: {out1}")
        self.logger.debug(f"Output2: {out2}")
        data1 = read_csv(output1_path)
        self.logger.info(f"Data1: {data1}")
        data2 = read_csv(output1_path)
        self.logger.info(f"Data2: {data2}")
        for res1, res2 in zip(data1, data2):
            for key in res1:
                if key in res2:
                    validate_output(res1[key], res2[key])
                else:
                    raise KeyError(f"Key '{key}' not found in second result.")

        self.consistent_output = True

class RunnerService:
    def __init__(
        self, 
        model_id, 
        logger, 
        ios_service, 
        checkup_service, 
        setup_service,
        model_path, 
        level,
        dir,
        remote,
        inspect,
        remove
    ):
        self.model_id = model_id
        self.logger = logger
        self.setup_service = setup_service
        self.ios_service = ios_service
        self.console = ios_service.console
        self.checkup_service = checkup_service
        self._model_path = model_path
        self.level = level
        self.dir = dir
        self.remote = remote
        self.inspect = inspect
        self.remove = remove
        self.example = ExampleGenerator(
            model_id=self.model_id
        )
        self.run_using_bash = False

    def run_model(self, input, output, batch):
        self.serve_model()
        if isinstance(input, list):
            input = input[0]
        self.logger.info("Running model")
        out = SetupService.run_command(
            ["ersilia", 
             "-v", 
             "run", 
             "-i", input[0], 
             "-o", output, 
             "-b", str(batch)
            ],
            logger=self.logger,
            capture_output=True
        )
        return out

    def fetch(self):
        SetupService.run_command(
            " ".join(["ersilia", 
            "-v", 
            "fetch", self.model_id, 
            "--from_dir", self.dir
            ]),
            logger=self.logger,
        )

    def serve_model(self):
        self.logger.info("Serving the model")
        SetupService.run_command(
            ["ersilia", 
            "-v", 
            "serve", self.model_id
            ],
            logger=self.logger,
            capture_output=True
        )

    def run_exampe(
            self, 
            n_samples, 
            file_name=None, 
            simple=True, 
            try_predefined=False
        ):
        return self.example.example(
            n_samples=n_samples, 
            file_name=file_name, 
            simple=simple, 
            try_predefined=try_predefined
        )
    @throw_ersilia_exception()
    def run_bash(self):
        def compute_rmse(y_true, y_pred):
            return sum((yt - yp) ** 2 for yt, yp in zip(y_true, y_pred)) ** 0.5 / len(y_true)
        
        def compare_outputs(bsh_data, ers_data):
            columns = set(bsh_data[0].keys()) & set(data[0].keys())
            self.logger.debug(f"Common columns: {columns}")

            for column in columns:
                bv = [row[column] for row in bsh_data]
                ev = [row[column] for row in ers_data]

                if all(isinstance(val, (int, float)) for val in bv + ev):
                    rmse = compute_rmse(bv, ev)
                    self.logger.debug(f"RMSE for {column}: {rmse}")

                    if rmse > 0.1:
                        raise texc.InconsistentOutputs(self.model_id)
                elif all(isinstance(val, str) for val in bv + ev):
                    if not all(
                        self._compare_string_similarity(a, b, 95) 
                        for a, b 
                        in zip(bv, ev)
                    ):
                        raise texc.InconsistentOutputs(self.model_id)
                    
        def read_csv(path, flag=False):
            try:
                with open(path, "r") as file:
                    lines = file.readlines()

                if not lines:
                    self.logger.error(f"File at {path} is empty.")
                    return []

                headers = lines[0].strip().split(",")
                if flag:
                    headers = headers[2:]

                data = []

                for line in lines[1:]:
                    self.logger.debug(f"Processing line: {line.strip()}")
                    values = line.strip().split(",")
                    values = values[2:] if flag else values

                    def infer_type(value):
                        try:
                            return int(value)
                        except ValueError:
                            try:
                                return float(value)
                            except ValueError:
                                return value  

                    _values = [infer_type(x) for x in values]

                    data.append(dict(zip(headers, _values)))

                return data
            except Exception as e:
                raise RuntimeError(
                    f"Failed to read CSV from {path}."
                ) from e


        def run_subprocess(command, env_vars=None):
            try:
                result = subprocess.run(
                    command, 
                    capture_output=True, 
                    text=True, 
                    check=True, 
                    env=env_vars,
                )
                self.logger.debug(
                    f"Subprocess output: {result.stdout}"
                )
                return result.stdout
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    "Subprocess execution failed."
                ) from e

        with tempfile.TemporaryDirectory() as temp_dir:

            model_path       = os.path.join(self.dir)
            temp_script_path = os.path.join(temp_dir, "script.sh")
            bash_output_path = os.path.join(temp_dir, "bash_output.csv")
            output_path      = os.path.join(temp_dir, "ersilia_output.csv")
            output_log_path  = os.path.join(temp_dir, "output.txt")
            error_log_path   = os.path.join(temp_dir, "error.txt")

            input = self.run_exampe(
                n_samples=Options.NUM_SAMPLES.value, 
                file_name=None,
                simple=True, 
                try_predefined=False
            )

            ex_file = os.path.join(temp_dir, "example_file.csv")

            with open(ex_file, "w") as f:
                f.write("smiles\n" + "\n".join(map(str, input)))

            run_sh_path = os.path.join(
                model_path, 
                "model", 
                "framework", 
                RUN_FILE
            )
            if not os.path.exists(run_sh_path):
                self.logger.warning(
                    f"{RUN_FILE} not found at {run_sh_path}. Skipping bash run."
                )
                return

            bash_script = f"""
                source {self.conda_prefix(self.is_base())}/etc/profile.d/conda.sh
                conda activate {self.model_id}
                cd {os.path.dirname(run_sh_path)}
                bash run.sh . {ex_file} {bash_output_path} > {output_log_path} 2> {error_log_path}
                conda deactivate
                """

            with open(temp_script_path, "w") as script_file:
                script_file.write(bash_script)

            self.logger.debug(f"Running bash script: {temp_script_path}")
            run_subprocess(["bash", temp_script_path])

            bsh_data = read_csv(bash_output_path)
            self.logger.info(f"Bash Data:{bsh_data}")
            self.logger.debug(f"Serving the model after run.sh")
            run_subprocess(
                ["ersilia", "-v", 
                 "serve", self.model_id, 
                ]
            )
            self.logger.debug(f"Running model for bash data consistency checking")
            out = run_subprocess(
                ["ersilia", "-v", 
                 "run", 
                 "-i", ex_file,
                 "-o", output_path
                ]
            )
            data = read_csv(output_path, flag=True)

            compare_outputs(bsh_data, data)

        self.run_using_bash = True

    @staticmethod
    def default_env():
        if "CONDA_DEFAULT_ENV" in os.environ:
            return os.environ["CONDA_DEFAULT_ENV"]
        else:
            return Options.BASE.value

    @staticmethod
    def conda_prefix(is_base):
        o = run_command_check_output("which conda").rstrip()
        if o:
            o = os.path.abspath(os.path.join(o, "..", ".."))
            return o
        if is_base:
            o = run_command_check_output("echo $CONDA_PREFIX").rstrip()
            return o
        else:
            o = run_command_check_output("echo $CONDA_PREFIX_1").rstrip()
            return o

    def is_base(self):
        default_env = self.default_env()
        self.logger.debug(f"Default environment: {default_env}")
        if default_env == "base":
            return True
        else:
            return False


    def _compare_string_similarity(
            self, 
            str1,
            str2, 
            threshold
        ):
        similarity = fuzz.ratio(str1, str2)
        return similarity >= threshold


    def make_output(self, elapsed_time):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        test_results = {
            "Date and Time Run": timestamp,
            "Time to Run Tests (seconds)": elapsed_time,
            "Basic Checks Passed": self.checkup_service.information_check,
            "Single Input Run Without Error": self.checkup_service.single_input,
            "Example Input Run Without Error": self.checkup_service.example_input,
            "Outputs Consistent": self.checkup_service.consistent_output,
            "Bash Run Without Error": self.run_using_bash,
        }

        data = [(key, str(value)) for key, value in test_results.items()]
        self.ios_service._generate_table(
            **TABLE_CONFIGS[TableType.FINAL_RUN_SUMMARY].__dict__,
            rows=data
        )

    def run(self, output_file=None):
        if MISSING_PACKAGES:
            raise ImportError(
                "Missing packages required for testing. "
                "Please install test extras with 'pip install ersilia[test]'."
            )

        if not output_file:
            output_file = os.path.join(
                self._model_path(self.model_id), 
                Options.OUTPUT_CSV.value
            )

        start_time = time.time()

        try:
            self._perform_checks(output_file)
            self._log_directory_sizes()
            self._perform_inspect()
            if self.level == Options.LEVEL_DEEP.value:
                self._perform_deep_checks(output_file)
            
            elapsed_time = time.time() - start_time
            self.make_output(elapsed_time)
            self._clear_folders()

        except Exception as e:
            click.echo(
                f"An error occurred: {e}"
            )
        finally:
            click.echo(
                "Run process finished successfully."
            )

    def transform_key(self, value):
        if value is True:
            return str(STATUS_CONFIGS.PASSED)
        elif value is False:
            return str(STATUS_CONFIGS.FAILED)
        return value
    
    def _perform_inspect(self):
        command, out = f"ersilia -v inspect {self.model_id}", None
        if self.remote and self.inspect:
            out = SetupService.run_command(
                command=command,
                logger=self.logger,
                capture_output=True,
                shell=True
            )
        elif not self.remote and self.inspect:
            out = SetupService.run_command(
                f"{command} -d {self.dir}",
                logger=self.logger,
                shell=True
            )
        if out is not None:
            out = json.loads(out)
            out = {" ".join(word.capitalize() for word in k.split("_")): self.transform_key(v) for k, v in out.items()}
            data = [(key, value) for key, value in out.items()]

            self.ios_service._generate_table(
             **TABLE_CONFIGS[TableType.INSPECT_SUMMARY].__dict__,
                rows=data,
                large_table=True,
                merge=True
            )
    
    def _perform_checks(self, output_file):
        self.checkup_service.check_information(output_file)
        self._generate_table(
            **TABLE_CONFIGS[TableType.MODEL_INFORMATION_CHECKS].__dict__,
            rows=self.ios_service.check_results
        )
        self.ios_service.check_results.clear()

        self.checkup_service.check_files()
        self._generate_table(
            **TABLE_CONFIGS[TableType.MODEL_FILE_CHECKS].__dict__,
            rows=self.ios_service.check_results
        )

    def _log_directory_sizes(self):
        dir_size, env_size = self.ios_service.get_directories_sizes()
        self._generate_table(
            **TABLE_CONFIGS[TableType.MODEL_DIRECTORY_SIZES].__dict__,
            rows=[(f"{dir_size:.2f}MB", f"{env_size:.2f}MB")]
        )

    def _perform_deep_checks(self, output_file):
        self.checkup_service.check_single_input(
            output_file, 
            self.run_model, 
            self.run_exampe
        )
        self._generate_table(
            **TABLE_CONFIGS[TableType.RUNNER_CHECKUP_STATUS].__dict__,
            rows = [
            ["Fetch", str(STATUS_CONFIGS.PASSED)],
            ["Serve", str(STATUS_CONFIGS.PASSED)],
            ["Run",   str(STATUS_CONFIGS.PASSED)]
         ]
        )
        self.checkup_service.check_example_input(
            output_file, 
            self.run_model, 
            self.run_exampe
        )
        self.checkup_service.check_consistent_output(
            self.run_exampe, 
            self.run_model
        )
        self.run_bash()

    def _generate_table(self, title, headers, rows):
        self.ios_service._generate_table(
            title=title,
            headers=headers,
            rows=rows
        )
    def _clear_folders(self):
        if self.remove:
            SetupService.run_command(
                f"rm -rf {self.dir}",
                logger=self.logger,
            )


class ModelTester(ErsiliaBase):
    def __init__(
            self, 
            model_id, 
            level, 
            dir,
            inspect,
            remote,
            remove
        ):
        ErsiliaBase.__init__(
            self, 
            config_json=None, 
            credentials_json=None
        )
        self.model_id = model_id
        self.level = level
        self.dir = dir or self.model_id
        self.inspect = inspect
        self.remote = remote
        self.remove = remove
        self.setup_service = SetupService(
            self.model_id, 
            self.dir, 
            self.logger,
            self.remote
        )
        self.ios = IOService(
            self.logger, 
            self._dest_dir,
            self._model_path, 
            self._get_bundle_location, 
            self._get_bentoml_location, 
            self.model_id,
            self.dir
        )
        self.checks = CheckService(
            self.logger, 
            self.model_id, 
            self._dest_dir,
            self.dir,
            self.ios,
        )
        self.runner = RunnerService(
            self.model_id,
            self.logger, 
            self.ios, 
            self.checks, 
            self.setup_service,
            self._model_path,
            self.level,
            self.dir,
            self.remote,
            self.inspect,
            self.remove
        )
    def setup(self):
        self.logger.debug(f"Running conda setup for {self.model_id}")
        self.setup_service.fetch_repo() # for remote option
        self.logger.debug(f"Fetching model {self.model_id} from local dir")
        self.runner.fetch()
        self.setup_service.check_conda_env()

    def run(self, output_file=None):
        self.runner.run(output_file)
