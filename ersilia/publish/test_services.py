# TODO adapt to input-type agnostic. For now, it works only with Compound input types.
import json
import os
import csv
import subprocess
import tempfile
import asyncio
import time
import click
import types
from collections import defaultdict
from datetime import datetime

from ersilia.utils.conda import SimpleConda
from .. import ErsiliaModel, throw_ersilia_exception
from ..hub.fetch.fetch import ModelFetcher
from ..cli import echo
from ..core.session import Session
from ..default import EOS, INFORMATION_FILE
from ..io.input import ExampleGenerator
from ..utils.exceptions_utils import test_exceptions as texc
from ..utils.terminal import run_command_check_output

# Check if we have the required imports in the environment
MISSING_PACKAGES = False
try:
    from scipy.stats import spearmanr
    from fuzzywuzzy import fuzz
    from rich.console import Console
    from rich.table import Table
    from rich.box import SIMPLE
except ImportError:
    MISSING_PACKAGES = True

RUN_FILE = "run.sh"
DATA_FILE = "data.csv"
NUM_SAMPLES = 5
BOLD = "\033[1m"
RESET = "\033[0m"
BASE = "base"

TEST_MESSAGES = {
        "pkg_err": "Missing packages required for testing. \
        Please install test extras with 'pip install ersilia[test]'."
}

class IOService:
    def __init__(
        self, logger,
        dest_dir,
        model_path, 
        bundle_path, 
        bentoml_path, 
        model_id,
        env,
        type
    ):
        self.logger = logger
        self.model_id = model_id
        self.env = env if env is not None else BASE
        self.model_size = 0
        self.console = Console()
        self.check_results = []
        self.type = type
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
            self.check_results.append((check_name, "[green]✔[/green]"))
            return True
        except Exception as e:
            self.logger.error(f"Check '{check_name,additional_info}' failed: {e}")
            self.check_results.append((check_name, "[red]✖[/red]"))
            return False
        
    def _generate_table(self, title, headers, rows):
        table = Table(
            title=title, 
            border_style="#FFC0CB"
        )
        for header in headers:
            table.add_column(
                header, 
                justify="center", 
                no_wrap=True,
                width=20
            )

        for row in rows:
            table.add_row(*[str(cell) for cell in row])
        # render
        self.console.print(table)

    def _get_file_requirements(self):

        if self.type.lower() == "bentoml":
            return [
                "src/service.py",
                "pack.py",
                "Dockerfile",
                "metadata.json",
                "README.md",
                "model/framework/run.sh",
                "LICENSE",
            ]
        elif self.type.lower() == "ersilia":
            return [
                "install.yml",
                "metadata.json",
                "metadata.yml",
                "model/framework/example/input.csv",
                "model/framework/example/output.csv",
                "README.md",
                "model/framework/run.sh",
                "LICENSE",
            ]
        else:
            raise ValueError(f"Unsupported model type: {self.type}")


    def _get_environment_location(self):
        conda = SimpleConda()
        python_path = conda.get_python_path_env(environment=self.env)
        return os.path.dirname(os.path.dirname(python_path))
    
    def read_information(self):
        file = os.path.join(
            self._dest_dir, 
            self.model_id, 
            INFORMATION_FILE
        )
        self.logger.info(f"Dest: {self._dest_dir}")
        if not os.path.exists(file):
            raise FileNotFoundError(f"Information file does not exist for model {self.model_id}")
        with open(file, "r") as f:
            return json.load(f)

    def set_model_size(self, directory):
        return sum(
            os.path.getsize(os.path.join(dirpath, filename))
            for dirpath, _, filenames in os.walk(directory)
            for filename in filenames
        )

    def print_output(self, result, output):
        def write_output(data):
            if output is not None:
                with open(output.name, "w") as file:
                    json.dump(data, file)
            else:
                print(json.dumps(data, indent=4))

        if isinstance(result, types.GeneratorType):
            for r in result:
                write_output(r if r is not None else "Something went wrong")
        else:
            print(result)

    def _get_environment_location(self):
        conda = SimpleConda()
        python_path = conda.get_python_path_env(environment=BASE)
        return os.path.dirname(os.path.dirname(python_path))


    @throw_ersilia_exception()
    def get_directories_sizes(self):
        self.logger.debug(BOLD + "Calculating model size... " + RESET)

        def format_size(size_in_bytes):
            if size_in_bytes < 1024:
                return f"{size_in_bytes}B"  # Bytes
            size_kb = size_in_bytes / 1024
            if size_kb < 1024:
                return f"{size_kb:.2f}KB"   # Kilobytes
            size_mb = size_kb / 1024
            if size_mb < 1024:
                return f"{size_mb:.2f}MB"   # Megabytes
            size_gb = size_mb / 1024
            return f"{size_gb:.2f}GB"       # Gigabytes

        def get_directory_size(directory, include_symlinks=True):
            if not directory or not os.path.exists(directory):
                return 0, defaultdict(int), defaultdict(int)

            _types, _sizes, total_sz = defaultdict(int), defaultdict(int), 0

            for dirpath, _, filenames in os.walk(directory):
                for filename in filenames:
                    filepath = os.path.join(dirpath, filename)

                    if include_symlinks or not os.path.islink(filepath):
                        try:
                            size = os.path.getsize(filepath)
                            total_sz += size
                            file_extension = os.path.splitext(filename)[1]
                            _types[file_extension] += 1
                            _sizes[file_extension] += size
                        except OSError as e:
                            self.logger.warning(f"Failed to access file {filepath}: {e}")

            return total_sz, _types, _sizes

        directories = {
            "dest_dir": (self._model_path(model_id=self.model_id), False),
            "bundle_dir": (self._bundle_path(model_id=self.model_id), False),
            "bentoml_dir": (self._bentoml_path(model_id=self.model_id), False),
            "env_dir": (self._get_environment_location(), True),
        }

        total_size = 0
        directory_sizes = {}

        for label, (directory, include_symlinks) in directories.items():
            size, file_types, _ = get_directory_size(
                directory=directory, 
                include_symlinks=include_symlinks
            )
            formatted_size = format_size(size)
            directory_sizes[label] = formatted_size
            total_size += size

        total_size_formatted = format_size(total_size)

        self.logger.debug(BOLD + "Model Size Breakdown:" + RESET)
        self.logger.debug(f"Total: {total_size_formatted}")

        self.logger.debug("Sizes of directories:")
        for label, size in directory_sizes.items():
            self.logger.debug(f"{label}: {size}")

        self.model_size = total_size_formatted
        return directory_sizes


class CheckService:
    def __init__(
        self, 
        logger, 
        model_id, 
        dest_dir, 
        dir, 
        ios,
        type
        ):
        self.logger = logger
        self.model_id = model_id
        self._dest_dir = dest_dir
        self.dir = dir
        self.type = type
        self._run_check = ios._run_check
        self._generate_table = ios._generate_table
        self._get_file_requirements = ios._get_file_requirements
        self.console = ios.console
        self.check_results = ios.check_results
        self.information_check = False
        self.single_input = False
        self.example_input = False
        self.consistent_output = False

    def _check_file_existence(self, file_path):
    
        if file_path == "metadata.json" and "ersilia" in self.type.lower():
            json_exists = os.path.exists(os.path.join(self.dir, "metadata.json"))
            yaml_exists = os.path.exists(os.path.join(self.dir, "metadata.yml"))
            if not (json_exists or yaml_exists):
                raise FileNotFoundError(f"Neither 'metadata.json' nor 'metadata.yml' found.")
        else:
            if not os.path.exists(os.path.join(self.dir, file_path)):
                raise FileNotFoundError(f"File '{file_path}' does not exist.")

    def check_files(self):
     
        self.logger.debug(f"Checking file requirements for {self.type} model.")
        fr = self._get_file_requirements()

        for fp in fr:
            self.logger.debug(f"Checking file: {fp}")
            self._run_check(
                self._check_file_existence, 
                None,
                f"File: {fp}", 
                fp
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
        valid_tasks = {
            "Classification",
            "Regression",
            "Generative",
            "Representation",
            "Similarity",
            "Clustering",
            "Dimensionality reduction",
        }

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

        invalid_tasks = [task for task in tasks if task not in valid_tasks]
        if invalid_tasks:
            raise texc.InvalidEntry(
                "Task", message=f"Invalid tasks: {', '.join(invalid_tasks)}"
            )

        self.logger.debug("All tasks are valid.")



    def _check_model_input(self, data):
        self.logger.debug("Checking model input...")
        valid_inputs = [{"Compound"}, {"Protein"}, {"Text"}]
        if set(data["card"]["Input"]) not in valid_inputs:
            raise texc.InvalidEntry("Input")

    def _check_model_input_shape(self, data):
        self.logger.debug("Checking model input shape...")
        valid_input_shapes = {
            "Single", 
            "Pair", 
            "List", 
            "Pair of Lists", 
            "List of Lists"
        }
        if data["card"]["Input Shape"] not in valid_input_shapes:
            raise texc.InvalidEntry("Input Shape")


    def _check_model_output(self, data):
        self.logger.debug("Checking model output...")
        valid_outputs = {
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

        raw_outputs = data.get("card", {}).get("Output", "")
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
            if output not in valid_outputs
        ]
        if invalid_outputs:
            raise texc.InvalidEntry(
                "Output", 
                message=f"Invalid outputs: {', '.join(invalid_outputs)}"
            )

        self.logger.debug("All outputs are valid.")


    def _check_model_output_type(self, data):
        self.logger.debug("Checking model output type...")
        valid_output_types = [{"String"}, {"Float"}, {"Integer"}]
        if set(data["card"]["Output Type"]) not in valid_output_types:
            raise texc.InvalidEntry("Output Type")

    def _check_model_output_shape(self, data):
        self.logger.debug("Checking model output shape...")
        valid_output_shapes = {
            "Single", 
            "List", 
            "Flexible List", 
            "Matrix", 
            "Serializable Object"
        }
        if data["card"]["Output Shape"] not in valid_output_shapes:
            raise texc.InvalidEntry("Output Shape")



    @throw_ersilia_exception()
    def check_information(self, output):
        self.logger.debug("Checking that model information is correct")
        self.logger.debug(
            BOLD + f"Beginning checks for {self.model_id} model information:" + RESET
        )
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
    def check_single_input(self, output, fetch_and_serve, run_model):

        input_smiles = "COc1ccc2c(NC(=O)Nc3cccc(C(F)(F)F)n3)ccnc2c1"

        self.logger.debug(BOLD + "Testing model on single smiles input...\n" + RESET)
        fetch_and_serve()
        result = run_model(
            input=input_smiles, 
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
            n_samples=NUM_SAMPLES, 
            file_name=None, 
            simple=True, 
            try_predefined=False
        )

        self.logger.debug(
            BOLD + "\nTesting model on input of 5 smiles given by 'example' command...\n" + RESET
        )
        self.logger.debug(f"This is the input: {input_samples}")

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
                if not os.path.exists(file_path):
                    raise FileNotFoundError(f"File not found: {file_path}")
                with open(file_path, mode='r') as csv_file:
                    reader = csv.DictReader(csv_file)
                    return [row for row in reader]
                
        def is_json(data):
            try:
                if isinstance(data, str):
                    json.loads(data)
                elif isinstance(data, list):
                    json.dumps(data) 
                return True
            except (ValueError, TypeError):
                return False
            


        self.logger.debug(BOLD + "\nConfirming model produces consistent output..." + RESET)

        input_samples = run_example(
            n_samples=NUM_SAMPLES, 
            file_name=None, 
            simple=True, 
            try_predefined=False
        )
        result1 = run_model(
            input=input_samples, 
            output="output1.csv", 
            batch=100
        )
        result2 = run_model(
            input=input_samples, 
            output="output2.csv", 
            batch=100
        )

        if is_json(result1) and is_json(result2):
            data1, data2 = result1, result2
        else:
            data1 = read_csv("output1.csv")
            data2 = read_csv("output2.csv")

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
        model_path, 
        env,
        type,
        level,
        dir
    ):
        self.model_id = model_id
        self.logger = logger
        self.ios_service = ios_service
        self.console = ios_service.console
        self.checkup_service = checkup_service
        self._model_path = model_path
        self.env = env
        self.type = type
        self.level = level
        self.dir = dir
        self._output_type = self.ios_service.read_information()["card"]["Output Type"]
        session = Session(config_json=None)
        service_class = session.current_service_class()
        self.model = ErsiliaModel(
                self.model_id, 
                service_class=service_class
        )
        self.example = ExampleGenerator(model_id=self.model_id)
        self.fetcher = ModelFetcher(repo_path=self.dir)
        self.run_using_bash = False

    def run_model(self, input, output, batch):
        return self.model.run(
            input=input,
            output=output,
            batch_size=batch
    )
    def fetch_and_serve(self):
        asyncio.run(self.fetcher.fetch(self.model_id
        ))
        self.model.serve()
    def serve_model(self, input, output, batch):
        return self.model.run(
            input=input,
            output=output,
            batch_size=batch
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
                    if not all(self._compare_string_similarity(a, b, 95) for a, b in zip(bv, ev)):
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

                    try:
                        if self._output_type == ["Float"]:
                            _values = [float(x) for x in values]
                        elif self._output_type == ["Integer"]:
                            _values = [int(x) for x in values]
                    except ValueError as e:
                        self.logger.warning(f"Value conversion error: {e}")
                    data.append(dict(zip(headers, _values)))
                
                return data
            except Exception as e:
                raise RuntimeError(f"Failed to read CSV from {path}.") from e

        def run_subprocess(command, env_vars=None):
            try:
                result = subprocess.run(
                    command, 
                    capture_output=True, 
                    text=True, 
                    check=True, 
                    env=env_vars
                )
                self.logger.debug(f"Subprocess output: {result.stdout}")
                return result.stdout
            except subprocess.CalledProcessError as e:
                raise RuntimeError("Subprocess execution failed.") from e

        with tempfile.TemporaryDirectory() as temp_dir:

            model_path       = os.path.join(self.dir)
            temp_script_path = os.path.join(temp_dir, "script.sh")
            bash_output_path = os.path.join(temp_dir, "bash_output.csv")
            output_path      = os.path.join(temp_dir, "ersilia_output.csv")
            output_log_path  = os.path.join(temp_dir, "output.txt")
            error_log_path   = os.path.join(temp_dir, "error.txt")

            input = self.run_exampe(
                n_samples=NUM_SAMPLES, 
                file_name=None,
                simple=True, 
                try_predefined=False
            )

            ex_file = os.path.join(temp_dir, "example_file.csv")

            with open(ex_file, "w") as example_file:
                example_file.write("smiles\n" + "\n".join(map(str, input)))

            run_sh_path = os.path.join(
                model_path, 
                "model", 
                "framework", 
                "run.sh"
            )
            if not os.path.exists(run_sh_path):
                self.logger.warning(f"run.sh not found at {run_sh_path}. Skipping bash run.")
                return

            bash_script = f"""
                source {self.conda_prefix(self.is_base())}/etc/profile.d/conda.sh
                conda activate {BASE}
                cd {os.path.dirname(run_sh_path)}
                bash run.sh . {ex_file} {bash_output_path} > {output_log_path} 2> {error_log_path}
                conda deactivate
            """

            with open(temp_script_path, "w") as script_file:
                script_file.write(bash_script)

            self.logger.debug(f"Running bash script: {temp_script_path}")
            run_subprocess(["bash", temp_script_path])

            bsh_data = read_csv(bash_output_path)
   
            self.run_model(
                ex_file, 
                output_path, 
                100
            )
            data = read_csv(output_path, flag=True)

            compare_outputs(bsh_data, data)

        self.run_using_bash = True

    @staticmethod
    def default_env():
        if "CONDA_DEFAULT_ENV" in os.environ:
            return os.environ["CONDA_DEFAULT_ENV"]
        else:
            return BASE

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


    def make_output(self, time, model_size):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.information_check = self.checkup_service.information_check
        self.single_input = self.checkup_service.single_input
        self.example_input = self.checkup_service.example_input
        self.consistent_output = self.checkup_service.consistent_output
        data = [
            ("Date and Time Run", timestamp),
            ("Model Size (MB)", model_size),
            ("Time to Run Tests (seconds)", time),
            ("Basic Checks Passed", self.information_check),
            ("Single Input Run Without Error", self.single_input),
            ("Example Input Run Without Error", self.example_input),
            ("Outputs Consistent", self.consistent_output),
            ("Bash Run Without Error", self.run_using_bash),
        ]
        
        headers = ["Check Type", "Status"]
        
        self.ios_service._generate_table("Test Run Summary", headers, data)

    def run(self, output_file=None):
        if MISSING_PACKAGES:
            raise ImportError(TEST_MESSAGES["pkg_err"])

        if not output_file:
            output_file = os.path.join(
                self._model_path(self.model_id), "result.csv"
            )

        st = time.time()
        try:
            self.checkup_service.check_information(output_file)
            self.ios_service._generate_table(
                title="Model Information Checks",
                headers=["Check", "Status"],
                rows=self.ios_service.check_results
            )
            self.ios_service.check_results.clear()
            self.checkup_service.check_files()
            self.ios_service._generate_table(
                title="Model Information Checks",
                headers=["Check", "Status"],
                rows=self.ios_service.check_results
            )
            ds = self.ios_service.get_directories_sizes()
            self.ios_service._generate_table(
                title="Model Directory Sizes",
                headers=["Dest dir", "Env Dir"],
                rows=[(ds["dest_dir"], ds["env_dir"])]
            )
            if self.level == "deep":
                self.checkup_service.check_single_input(
                    output_file, 
                    self.fetch_and_serve,
                    self.run_model
                )
                self.ios_service._generate_table(
                    title="Runner Checkup Status",
                    headers=["Runner", "Status"],
                    rows=[["Fetch", "[green]✔[/green]"], ["Serve", "[green]✔[/green]"], ["Run", "[green]✔[/green]"]]
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
                model_size = self.ios_service.model_size
                self.run_bash()
            et = time.time()
            elapsed = et - st
            self.make_output(elapsed, model_size)

        except Exception as e:
            click.echo(f"❌ An error occurred: {e}")
        finally:
            click.echo("🏁 Run process finished.")