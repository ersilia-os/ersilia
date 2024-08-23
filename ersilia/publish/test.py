# TODO adapt to input-type agnostic. For now, it works only with Compound input types.
import json
import os
import subprocess
import tempfile
import time
import click
import types
from collections import defaultdict
from datetime import datetime

import numpy as np

from ersilia.utils.conda import SimpleConda
from .. import ErsiliaBase, ErsiliaModel, throw_ersilia_exception
from ..cli import echo
from ..core.session import Session
from ..default import EOS, INFORMATION_FILE
from ..io.input import ExampleGenerator
from ..utils.exceptions_utils import test_exceptions as texc
from ..utils.logging import make_temp_dir
from ..utils.terminal import run_command_check_output



# Check if we have the required imports in the environment
MISSING_PACKAGES = False
try:
    from scipy.stats import spearmanr
    from sklearn.metrics import mean_squared_error
    from fuzzywuzzy import fuzz
except ImportError:
    MISSING_PACKAGES = True

RUN_FILE = "run.sh"
DATA_FILE = "data.csv"
NUM_SAMPLES = 5
BOLD = "\033[1m"
RESET = "\033[0m"


class ModelTester(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.model_size = 0
        self.tmp_folder = make_temp_dir(prefix="ersilia-")
        self._info = self._read_information()
        self._input = self._info["card"]["Input"]
        self._output_type = self._info["card"]["Output Type"]
        self.RUN_FILE = "run.sh"
        self.information_check = False
        self.single_input = False
        self.example_input = False
        self.consistent_output = False
        self.run_using_bash = False

    def _read_information(self):
        json_file = os.path.join(self._dest_dir, self.model_id, INFORMATION_FILE)
        self.logger.debug("Reading model information from {0}".format(json_file))
        if not os.path.exists(json_file):
            raise texc.InformationFileNotExist(self.model_id)
        with open(json_file, "r") as f:
            data = json.load(f)
        return data

    """
    This function uses the fuzzy wuzzy package to compare the differences between outputs when 
    they're strings and not floats. The fuzz.ratio gives the percent of similarity between the two outputs.
    Example: two strings that are the exact same will return 100
    """

    def _compare_output_strings(self, output1, output2):
        if output1 is None and output2 is None:
            return 100
        else:
            return fuzz.ratio(output1, output2)

   
   
    """
    When the user specifies an output file, the file will show the user how big the model is. This function 
    calculates the size of the model to allow this. 
    """

    def _set_model_size(self, directory):
        for dirpath, dirnames, filenames in os.walk(directory):
            for filename in filenames:
                file_path = os.path.join(dirpath, filename)
                self.model_size += os.path.getsize(file_path)

    """
    This helper method was taken from the run.py file, and just prints the output for the user 
    """

    def _print_output(self, result, output):
        echo("Printing output...")

        if isinstance(result, types.GeneratorType):
            for r in result:
                if r is not None:
                    if output is not None:
                        with open(output.name, "w") as file:
                            json.dump(r, output.name)
                    else:
                        echo(json.dumps(r, indent=4))
                else:
                    if output is not None:
                        message = echo("Something went wrong", fg="red")
                        with open(output.name, "w") as file:
                            json.dump(message, output.name)
                    else:
                        echo("Something went wrong", fg="red")

        else:
            echo(result)

    """
    This helper method checks that the model ID is correct.
    """

    def _check_model_id(self, data):
        self.logger.debug("Checking model ID...")
        if data["card"]["Identifier"] != self.model_id:
            raise texc.WrongCardIdentifierError(self.model_id)

    """
    This helper method checks that the slug field is non-empty.
    """

    def _check_model_slug(self, data):
        self.logger.debug("Checking model slug...")
        if not data["card"]["Slug"]:
            raise texc.EmptyField("slug")

    """
    This helper method checks that the description field is non-empty.
    """

    def _check_model_description(self, data):
        self.logger.debug("Checking model description...")
        if not data["card"]["Description"]:
            raise texc.EmptyField("Description")

    """
    This helper method checks that the model task is one of the following valid entries:
        - Classification
        - Regression
        - Generative
        - Representation
        - Similarity
        - Clustering
        - Dimensionality reduction
    """

    def _check_model_task(self, data):
        self.logger.debug("Checking model task...")
        valid_tasks = [
            "Classification",
            "Regression",
            "Generative",
            "Representation",
            "Similarity",
            "Clustering",
            "Dimensionality reduction",
        ]
        sep = ", "
        tasks = []
        if sep in data["card"]["Task"]:
            tasks = data["card"]["Task"].split(sep)
        else:
            tasks = data["card"]["Task"]
        for task in tasks:
            if task not in valid_tasks:
                raise texc.InvalidEntry("Task")

    """
    This helper method checks that the input field is one of the following valid entries:
        - Compound
        - Protein
        - Text
    """

    def _check_model_input(self, data):
        self.logger.debug("Checking model input...")
        valid_inputs = [["Compound"], ["Protein"], ["Text"]]
        if data["card"]["Input"] not in valid_inputs:
            raise texc.InvalidEntry("Input")

    """
    This helper method checks that the input shape field is one of the following valid entries:
        - Single
        - Pair
        - List
        - Pair of Lists
        - List of Lists
    """

    def _check_model_input_shape(self, data):
        self.logger.debug("Checking model input shape...")
        valid_input_shapes = [
            "Single",
            "Pair",
            "List",
            "Pair of Lists",
            "List of Lists",
        ]
        if data["card"]["Input Shape"] not in valid_input_shapes:
            raise texc.InvalidEntry("Input Shape")

    """
    This helper method checks the the output is one of the following valid entries:
        - Boolean
        - Compound
        - Descriptor
        - Distance
        - Experimental value
        - Image
        - Other value
        - Probability
        - Protein
        - Score
        - Text
    """

    def _check_model_output(self, data):
        self.logger.debug("Checking model output...")
        valid_outputs = [
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
        ]
        sep = ", "
        outputs = []
        if sep in data["card"]["Output"]:
            outputs = data["card"]["Output"].split(sep)
        else:
            outputs = data["card"]["Output"]
        for output in outputs:
            if output not in valid_outputs:
                raise texc.InvalidEntry("Output")

    """
    This helper method checks that the output type is one of the following valid entries:
        - String
        - Float
        - Integer
    """

    def _check_model_output_type(self, data):
        self.logger.debug("Checking model output type...")
        valid_output_types = [["String"], ["Float"], ["Integer"]]
        if data["card"]["Output Type"] not in valid_output_types:
            raise texc.InvalidEntry("Output Type")

    """
    This helper method checks that the output shape is one of the following valid entries:
        - Single
        - List
        - Flexible List
        - Matrix
        - Serializable Object
    """

    def _check_model_output_shape(self, data):
        self.logger.debug("Checking model output shape...")
        valid_output_shapes = [
            "Single",
            "List",
            "Flexible List",
            "Matrix",
            "Serializable Object",
        ]
        if data["card"]["Output Shape"] not in valid_output_shapes:
            raise texc.InvalidEntry("Output Shape")

    """
    Check the model information to make sure it's correct. Performs the following checks:
    - Checks that model ID is correct
    - Checks that model slug is non-empty
    - Checks that model description is non-empty
    - Checks that the model task is valid
    - Checks that the model input, input shape is valid
    - Checks that the model output, output type, output shape is valid
    """

    @throw_ersilia_exception
    def check_information(self, output):


        self.logger.debug("Checking that model information is correct")
        self.logger.debug(
            BOLD
            + "Beginning checks for {0} model information:".format(self.model_id)
            + RESET
        )
        json_file = os.path.join(self._dest_dir, self.model_id, INFORMATION_FILE)
        with open(json_file, "r") as f:
            data = json.load(f)

        self._check_model_id(data)
        self._check_model_slug(data)
        self._check_model_description(data)
        self._check_model_task(data)
        self._check_model_input(data)
        self._check_model_input_shape(data)
        self._check_model_output(data)
        self._check_model_output_type(data)
        self._check_model_output_shape(data)
        click.echo(BOLD + "Test: Model information, SUCCESS! ✅\n" + RESET)

        if output is not None:
            self.information_check = True
        

    """
    Runs the model on a single smiles string and prints to the user if no output is specified.
    """

    @throw_ersilia_exception
    def check_single_input(self, output):
        session = Session(config_json=None)
        service_class = session.current_service_class()
        input = "COc1ccc2c(NC(=O)Nc3cccc(C(F)(F)F)n3)ccnc2c1"

        self.logger.debug(BOLD + "Testing model on single smiles input...\n" + RESET)
        mdl = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        result = mdl.run(input=input, output=output, batch_size=100)

        if output is not None:
            self.single_input = True
            click.echo(BOLD + "Test: Single SMILES input, SUCCESS! ✅\n" + RESET)
        else:
            self._print_output(result, output)

    """
    Generates an example input of 5 smiles using the 'example' command, and then tests the model on that input and prints it
    to the consol if no output file is specified by the user.
    """

    @throw_ersilia_exception
    def check_example_input(self, output):
        session = Session(config_json=None)
        service_class = session.current_service_class()
        eg = ExampleGenerator(model_id=self.model_id)
        input = eg.example(n_samples=NUM_SAMPLES, file_name=None, simple=True, try_predefined=False) 
        self.logger.debug(
            BOLD
            + "\nTesting model on input of 5 smiles given by 'example' command...\n"
            + RESET
        )
        self.logger.debug("This is the input: {0}".format(input))
        mdl = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        result = mdl.run(input=input, output=output, batch_size=100)
        if output is not None:
            self.example_input = True
            click.echo(BOLD + "Test: 5 SMILES input, SUCCESS! ✅\n")
        else:
            self._print_output(result, output)

    """
    Gets an example input of 5 smiles using the 'example' command, and then runs this same input on the 
    model twice. Then, it checks if the outputs are consistent or not and specifies that to the user. If 
    it is not consistent, an InconsistentOutput error is raised. Lastly, it makes sure that the number of 
    outputs equals the number of inputs.  
    """
    
    
    
    
    
    
    @throw_ersilia_exception
    def check_consistent_output(self):
        def compute_rmse(y_true, y_pred):
            squared_errors = [(yt - yp) ** 2 for yt, yp in zip(y_true, y_pred)]
            mse = sum(squared_errors) / len(squared_errors)
            return mse

    
        self.logger.debug(BOLD + "\nConfirming model produces consistent output..." + RESET)

        session = Session(config_json=None)
        service_class = session.current_service_class()

        eg = ExampleGenerator(model_id=self.model_id)
        input = eg.example(n_samples=NUM_SAMPLES, file_name=None, simple=True, try_predefined=False)
        input = eg.example(n_samples=NUM_SAMPLES, file_name=None, simple=True, try_predefined=False)

        mdl1 = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        mdl2 = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        mdl1 = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        mdl2 = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        result = mdl1.run(input=input, output=None, batch_size=100)
        result2 = mdl2.run(input=input, output=None, batch_size=100)

        zipped = list(zip(result, result2))

        for item1, item2 in zipped:
            output1 = item1["output"]
            output2 = item2["output"]

            keys1 = list(output1.keys())
            keys2 = list(output2.keys())

            for key1, key2 in zip(keys1, keys2):
                if not isinstance(output1[key1], type(output2[key2])):
                    for item1, item2 in zipped:
                        self.logger.debug(item1)
                        self.logger.debug(item2)
                        self.logger.debug("\n")
                        self.logger.debug(item1)
                        self.logger.debug(item2)
                        self.logger.debug("\n")
                    raise texc.InconsistentOutputTypes(self.model_id)

                if output1[key1] is None:
                    continue

                elif isinstance(output1[key1], (float, int)):
                    # Calculate RMSE
                    rmse = compute_rmse([output1[key1]], [output2[key2]])
                    self.logger.debug(f"RMSE for {key1}: {rmse}")
                    if rmse > 0.1:  # Adjust the threshold as needed
                        self.logger.debug(
                            BOLD
                            + "\nBash run and Ersilia run produce inconsistent results (Root Mean Square Error difference exceeds 10%)."
                            + RESET
                        )
                        raise texc.InconsistentOutputs(self.model_id)

                    # Calculate Spearman's correlation
                    rho, p_value = spearmanr([output1[key1]], [output2[key2]])
                    self.logger.debug(f"Spearman's correlation for {key1}: {rho}")
                    if rho < 0.5:  # Adjust the threshold as needed
                        self.logger.debug(
                            BOLD
                            + "\nBash run and Ersilia run produce inconsistent results (Spearman's correlation below threshold)."
                            + RESET
                        )
                        raise texc.InconsistentOutputs(self.model_id)

                    elif isinstance(output1[key1], list):
                        ls1 = output1[key1]
                        ls2 = output2[key2]

                        # Calculate rmse for lists
                        rmse = compute_rmse(ls1, ls2)
                        self.logger.debug(f"RMSE for {key1}: {rmse}")
                        if rmse > 0.1:  # Adjust the threshold as needed
                            self.logger.debug(
                                BOLD
                                + "\nBash run and Ersilia run produce inconsistent results (Root Mean Square Error exceeded threshold of 10% for list)."
                                + RESET
                            )
                            raise texc.InconsistentOutputs(self.model_id)

                        # Calculate Spearman's correlation for lists
                        rho, p_value = spearmanr(ls1, ls2)
                        self.logger.debug(f"Spearman's correlation for {key1}: {rho}")
                        if rho < 0.5:  # Adjust the threshold as needed
                            self.logger.debug(
                                BOLD
                                + "\nBash run and Ersilia run produce inconsistent results (Spearman's correlation below threshold for list)."
                                + RESET
                            )
                            raise texc.InconsistentOutputs(self.model_id)
                    # Check String Outputs
                    else: 
                        if self._compare_output_strings(output1[key1], output2[key2]) <= 95:
                            self.logger.debug(f"output1 value: {output1[key1]}")
                            self.logger.debug(f"output2 value: {output2[key2]}")
                            raise texc.InconsistentOutputs(self.model_id)
        self.consistent_output = True
        click.echo(
            BOLD + "Test: Output Consistency, SUCCESS! ✅\n" + RESET
        )

        self.logger.debug(
            BOLD + "\nConfirming there are same number of outputs as inputs..." + RESET
        )
        self.logger.debug(f"Number of inputs: {NUM_SAMPLES}")
        self.logger.debug(f"Number of outputs: {len(zipped)}")

        if NUM_SAMPLES != len(zipped):
            self.logger.debug("Number of inputs: {NUM_SAMPLES}")
            self.logger.debug("Number of outputs: {len(zipped)}")
            raise texc.MissingOutputs()
        else:
            click.echo(BOLD + "\nTest: Equal Inputs and Outputs, SUCCESS! ✅\n" + RESET)

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
        if default_env == "base":
            return True
        else:
            return False


    def _compare_string_similarity(self, str1, str2, similarity_threshold):
        similarity = fuzz.ratio(str1, str2)
        return similarity >= similarity_threshold

    
    @staticmethod
    def get_directory_size_without_symlinks(directory):
        if directory is None:
            return 0, {}, {}
        if not os.path.exists(directory):
            return 0, {}, {}
        
        total_size = 0
        file_types = defaultdict(int)
        file_sizes = defaultdict(int)
        
        for dirpath, dirnames, filenames in os.walk(directory):
            for filename in filenames:
                filepath = os.path.join(dirpath, filename)
                if not os.path.islink(filepath):
                    size = os.path.getsize(filepath)
                    total_size += size
                    file_extension = os.path.splitext(filename)[1]
                    file_types[file_extension] += 1
                    file_sizes[file_extension] += size
        
        return total_size, file_types, file_sizes


    @staticmethod
    def get_directory_size_with_symlinks(directory):
        if directory is None:
            return 0, {}, {}
        if not os.path.exists(directory):
            return 0, {}, {}
        
        total_size = 0
        file_types = defaultdict(int)
        file_sizes = defaultdict(int)
        
        for dirpath, dirnames, filenames in os.walk(directory):
            for filename in filenames:
                filepath = os.path.join(dirpath, filename)
                size = os.path.getsize(filepath)
                total_size += size
                file_extension = os.path.splitext(filename)[1]
                file_types[file_extension] += 1
                file_sizes[file_extension] += size
        return total_size, file_types, file_sizes
        
    # Get location of conda environment
    def _get_environment_location(self):
        conda = SimpleConda()
        python_path = conda.get_python_path_env(environment=self.model_id)
        env_dir = os.path.dirname(python_path).split("/")
        env_dir = "/".join(env_dir[:-1])
        return env_dir
    

    @throw_ersilia_exception
    def get_directories_sizes(self):
        self.logger.debug(BOLD + "Calculating model size... " + RESET)
        def log_file_analysis(size, file_types, file_sizes, label):
            self.logger.debug(f"Analyzing files in {label}:")
            self.logger.debug(f"File types & counts: {dict(file_types)}")
            self.logger.debug(f"Total size: {size} bytes")
        dest_dir = self._model_path(model_id=self.model_id)
        bundle_dir = self._get_bundle_location(model_id=self.model_id)
        bentoml_dir = self._get_bentoml_location(model_id=self.model_id)
        env_dir = self._get_environment_location()
        dest_size, dest_file_types, dest_file_sizes = self.get_directory_size_without_symlinks(dest_dir)
        bundle_size, bundle_file_types, bundle_file_sizes = self.get_directory_size_without_symlinks(bundle_dir)
        bentoml_size, bentoml_file_types, bentoml_file_sizes = self.get_directory_size_without_symlinks(bentoml_dir)
        env_size, env_file_types, env_file_sizes = self.get_directory_size_with_symlinks(env_dir)

        log_file_analysis(dest_size, dest_file_types, dest_file_sizes, "dest_dir")
        log_file_analysis(bundle_size, bundle_file_types, bundle_file_sizes, "bundle_dir")
        log_file_analysis(bentoml_size, bentoml_file_types, bentoml_file_sizes, "bentoml_dir")
        log_file_analysis(env_size, env_file_types, env_file_sizes, "env_dir")
        
        model_size = dest_size + bundle_size + bentoml_size + env_size
        self.model_size = model_size
        size_kb = model_size / 1024
        size_mb = size_kb / 1024
        size_gb = size_mb / 1024

        self.logger.debug(BOLD + "Model Size Breakdown:" + RESET)
        self.logger.debug(f"KB: {size_kb}")
        self.logger.debug(f"MB: {size_mb}")
        self.logger.debug(f"GB: {size_gb}")

        self.logger.debug("Sizes of directories:")
        self.logger.debug(f"dest_size: {dest_size} bytes")
        self.logger.debug(f"bundle_size: {bundle_size} bytes")
        self.logger.debug(f"bentoml_size: {bentoml_size} bytes")
        self.logger.debug(f"env_size: {env_size} bytes")
        return {
        "dest_size": dest_size,
        "bundle_size": bundle_size,
        "bentoml_size": bentoml_size,
        "env_size": env_size
    }
    
    
    @throw_ersilia_exception
    def run_bash(self):
        def updated_read_csv(self, file_path, ersilia_flag = False):
            data = []
            with open(file_path, "r") as file:
                lines = file.readlines()
                headers = lines[0].strip().split(",")
                if ersilia_flag:
                    headers = headers[2:]
                
                print("\n", "\n")
                
                for line in lines[1:]:
                    self.logger.debug(f"Processing line: {line}")
                    values = line.strip().split(",")
                    if ersilia_flag:
                        selected_values = values[2:]
                    else:
                        selected_values = values
                    self.logger.debug(f"Selected Values: {selected_values} and their type {self._output_type}")
                    
                    if self._output_type == ["Float"]:
                        selected_values = [float(x) for x in selected_values]
                        self.logger.debug(f"Converted to floats: {selected_values}")
                    elif self._output_type == ["Integer"]:
                        selected_values = [int(x) for x in selected_values]
                        self.logger.debug(f"Converted to integers: {selected_values}")
                    else:
                        self.logger.debug(f"Unknown type, keeping as strings: {selected_values}")
                    
                    row_data = dict(zip(headers, selected_values))
                    self.logger.debug(f"Appending row data: {row_data}")
                    data.append(row_data)
            
            return data
        
        with tempfile.TemporaryDirectory() as temp_dir:
            self.logger.debug(BOLD + "\nRunning the model bash script..." + RESET)  
            model_path =  os.path.join(EOS, "dest", self.model_id)

            # Create an example input
            eg = ExampleGenerator(model_id=self.model_id)
            input = eg.example(n_samples=NUM_SAMPLES, file_name=None, simple=True, try_predefined=False)
            # Read it into a temp file
            ex_file = os.path.abspath(os.path.join(temp_dir, "example_file.csv"))
            
            
            with open(ex_file, "w") as f:
                f.write("smiles")
                for item in input:
                    f.write(str(item) + "\n")

            run_sh_path = os.path.join(model_path, "model", "framework", "run.sh")
            run_sh_path = os.path.join(model_path, "model", "framework", "run.sh")
            # Halt this check if the run.sh file does not exist (e.g. eos3b5e)
            if not os.path.exists (run_sh_path):
                self.logger.debug(
                    BOLD + "\n ⛔️ Check halted. Either run.sh file does not exist, or model was not fetched via --from_github or --from_s3. ⛔️" + RESET
                )
                return

            # Navigate into the temporary directory
           
            subdirectory_path = os.path.join(model_path, "model", "framework")
            self.logger.debug(f"Changing directory to: {subdirectory_path}")
            os.chdir(subdirectory_path)
            try:
                run_path = os.path.abspath(subdirectory_path)
                tmp_script = os.path.abspath(os.path.join(temp_dir, "script.sh"))
                bash_output_path = os.path.abspath(os.path.join(temp_dir, "bash_output.csv"))
                output_log = os.path.abspath(os.path.join(temp_dir, "output.txt"))
                error_log = os.path.abspath(os.path.join(temp_dir, "error.txt"))
                bash_script = """
    source {0}/etc/profile.d/conda.sh     
    source {0}/etc/profile.d/conda.sh     
    conda activate {1}
    cd {2}
    bash run.sh . {3} {4} > {5} 2> {6}
    conda deactivate
    """.format(
                    self.conda_prefix(self.is_base()),
                    self.model_id,
                    run_path,
                    ex_file,
                    bash_output_path,
                    bash_output_path,
                    output_log,
                    error_log,
                )
                self.logger.debug(f"Script path: {tmp_script}")
                self.logger.debug(f"bash output path: {bash_output_path}")
                self.logger.debug(f"Output log path: {output_log}")
                self.logger.debug(f"Error log path: {error_log}")
                with open(tmp_script, "w") as f:
                    f.write(bash_script)

                

                self.logger.debug(BOLD + "\nExecuting 'bash run.sh'..." + RESET)
                try:
                    bash_result = subprocess.run(
                        ["bash", tmp_script], capture_output=True, text=True, check=True
                    )
                except subprocess.CalledProcessError as e:
                    self.logger.debug(f"STDOUT: {e.stdout}")
                    raise RuntimeError("Error encountered while running the bash script.") from e

                if os.path.exists(bash_output_path):
                    with open(bash_output_path, "r") as bash_output_file:
                        output_content = bash_output_file.read()
                        self.logger.debug(BOLD + "\nBash execution completed!" + RESET)
                        self.logger.debug("Captured Raw Bash Output:")
                        self.logger.debug(output_content)
                else:
                    click.echo(BOLD + "\nWARNING: Bash output file not found when reading the path:"+RESET)
                    self.logger.debug(bash_output_path)
                    click.echo(BOLD + "\n Ersilia and Bash comparison will raise an error" + RESET)
                    self.logger.debug(f"Generating bash script content for debugging: {tmp_script}\n")
             
                with open(error_log, "r") as error_file:
                    error_content = error_file.read()
                    self.logger.debug(BOLD + "\nCaptured Error:" + RESET)
                    if error_content == "":
                        click.echo("No errors on bash run found 😄 ✅\n")
                        self.run_using_bash = True # bash run was successful
                    else:
                        self.logger.debug(error_content)

            except Exception as e:
                raise RuntimeError(f"Error while activating the conda environment: {e}")

            self.logger.debug(BOLD + "\nExecuting ersilia run..."+RESET)
            ersilia_output_path = os.path.abspath(os.path.join(temp_dir, "ersilia_output.csv"))
            self.logger.debug(f"Ersilia output will be written to: {ersilia_output_path}")

            session = Session(config_json=None)
            service_class = session.current_service_class()
            mdl = ErsiliaModel(
                self.model_id, service_class=service_class, config_json=None
            )
            result = mdl.run(input=ex_file, output=ersilia_output_path, batch_size=100) 
           


            if os.path.exists(ersilia_output_path):
                with open(ersilia_output_path, "r") as ersilia_output_file:
                    output_content = ersilia_output_file.read()
                    click.echo("No errors on Ersilia run found 😄 ✅")
                    self.logger.debug("Captured Raw Ersilia Output:")
                    self.logger.debug(output_content)
            else:
                self.logger.debug(BOLD+f"Ersilia output file not found from the path: {ersilia_output_path}"+RESET)
            self.logger.debug("Processing ersilia csv output...")
            ersilia_run = updated_read_csv(self, ersilia_output_path, True)


            with open(bash_output_path, "r") as bash_output_file:
                    output_content = bash_output_file.read()
                    self.logger.debug("Captured Raw Bash Output:")
                    self.logger.debug(output_content)
            self.logger.debug("Processing raw bash output...: ")
            bash_run = updated_read_csv(self, bash_output_path, False)
            self.logger.debug(f"\nBash output:\n {bash_run}")
            self.logger.debug(f"\nErsilia output:\n {ersilia_run}")

            # Select common columns for comparison
            ersilia_columns = set()
            for row in ersilia_run:
                ersilia_columns.update(row.keys())
            self.logger.debug(f"\n Ersilia columns:  {ersilia_columns}")

            bash_columns = set()
            for row in bash_run:
                bash_columns.update(row.keys())
            self.logger.debug(f"\n Bash columns:  {bash_columns}")

            common_columns = ersilia_columns & bash_columns
            def compute_rmse(y_true, y_pred):
                squared_errors = [(yt - yp) ** 2 for yt, yp in zip(y_true, y_pred)]
                mse = sum(squared_errors) / len(squared_errors)
                return mse


 
            idx = 1
            self.logger.debug(BOLD + "\nComparing outputs from Ersilia and Bash runs..." + RESET)
            for column in common_columns:
                for i in range(len(ersilia_run)):
                    if isinstance(ersilia_run[i][column], (float, int)) and isinstance(
                        bash_run[i][column], (float, int)
                    ):
                        values1 = [row[column] for row in ersilia_run]
                        values2 = [row[column] for row in bash_run]
                        rmse = compute_rmse(values1, values2)
                        self.logger.debug(f"Root Mean Square Error for {column}: {rmse}")
                        if rmse > 0.10:  
                            self.logger.debug(
                                BOLD
                                + "\nBash run and Ersilia run produce inconsistent results (Root Mean Square Error difference exceeds 10%)."
                                + RESET
                            )
                            self.logger.debug(f"Values that raised error: {values1}, {values2}")
                            raise texc.InconsistentOutputs(self.model_id)
                        rho, p_value = spearmanr(values1, values2)
                        self.logger.debug(f"Spearman's correlation for {column}: {rho}\n")
                        if rho < 0.5:  # Adjust the threshold as needed
                            self.logger.debug(
                                BOLD
                                + "\nBash run and Ersilia run produce inconsistent results (Spearman's correlation below threshold)."
                                + RESET
                            )
                            raise texc.InconsistentOutputs(self.model_id)
                    # Both instances are strings
                    elif isinstance(ersilia_run[i][column], str) and isinstance(
                        bash_run[i][column], str
                    ):
                        if not all(
                            self._compare_string_similarity(a, b, 95)
                            for a, b in zip(ersilia_run[i][column], bash_run[i][column])
                        ):
                            self.logger.debug(
                                BOLD
                                + "\nBash run and Ersilia run produce inconsistent results."
                                + RESET
                            )
                            self.logger.debug(f"Error in the following column:  {column}")
                            self.logger.debug(ersilia_run[i][column])
                            self.logger.debug(bash_run[i][column])
                            raise texc.InconsistentOutputs(self.model_id)
                    elif isinstance(ersilia_run[i][column], bool) and isinstance(
                        ersilia_run[i][column], bool
                    ):
                        if not ersilia_run[i][column].equals(bash_run[i][column]):
                            self.logger.debug(
                                BOLD
                                + "\nBash run and Ersilia run produce inconsistent results."
                                + RESET
                            )
                            self.logger.debug("Error in the following column: ", column)
                            self.logger.debug(ersilia_run[i][column])
                            self.logger.debug(bash_run[i][column])
                            raise texc.InconsistentOutputs(self.model_id)

            click.echo(
                BOLD
                + "Test: Bash and Ersilia run comparison check, SUCCESS! ✅  Test Complete 🎉!"
                + RESET
            )

    """
    writes to the .json file all the basic information received from the test module:
    - size of the model
    - did the basic checks pass? True or False
    - time to run the model
    - did the single input run without error? True or False
    - did the run bash run without error? True or False
    - did the example input run without error? True or False 
    - are the outputs consistent? True or False 
    """

    def make_output(self, output, time):
        size_kb = self.model_size / 1024
        size_mb = size_kb / 1024
        size_gb = size_mb / 1024

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")


        data = {
            "date and time run": timestamp,  # Add date and time field
            "model size": {"KB": size_kb, "MB": size_mb, "GB": size_gb},
            "time to run tests (seconds)": time,
            "basic checks passed": self.information_check,
            "single input run without error": self.single_input,
            "example input run without error": self.example_input,
            "outputs consistent": self.consistent_output,
            "bash run without error": self.run_using_bash,
        }
        with open(output, "w") as json_file:
            json.dump(data, json_file, indent=4)



    def run(self, output_file):
        output_file = os.path.join(self._model_path(self.model_id), "TEST_MODULE_OUTPUT.csv") 
        start = time.time()
        self.check_information(output_file)
        self.check_single_input(output_file)
        self.check_example_input(output_file)
        self.check_consistent_output()
        self.get_directories_sizes()
        self.get_directories_sizes()
        self.run_bash()
        end = time.time()
        seconds_taken = end - start
        
        
        if output_file:
            self.make_output(output_file, seconds_taken)
            click.echo(f"📁 The output file is located at: {output_file}")
        else:
            click.echo("No output file specified! Skipping output file generation.")
