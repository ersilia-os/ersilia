import os
import json
import click
import tempfile
import types
import subprocess
from ..cli import echo

from ..io.input import ExampleGenerator
from .. import ErsiliaBase
from .. import throw_ersilia_exception
from .. import ErsiliaModel

from ..utils.exceptions_utils import test_exceptions as texc
from ..core.session import Session

from ..default import INFORMATION_FILE

RUN_FILE = "run.sh"
DATA_FILE = "data.csv"
OUTPUT_FILE = "output.csv"
FRAMEWORK_BASEDIR = "framework"
DIFFERENCE_THRESHOLD = 5    # outputs should be within this percent threshold to be considered consistent
NUM_SAMPLES = 5
BOLD = '\033[1m'
RESET = '\033[0m'
RED = '\033[31m'

class ModelTester(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        self._info = self._read_information()
        self._input = self._info["card"]["Input"]
        self._prepare_input_files()
        self.RUN_FILE = "run.sh"

    def _read_information(self):
        json_file = os.path.join(self._dest_dir, self.model_id, INFORMATION_FILE)
        self.logger.debug("Reading model information from {0}".format(json_file))
        if not os.path.exists(json_file): 
            raise texc.InformationFileNotExist(self.model_id)
        with open(json_file, "r") as f:
            data = json.load(f)
        return data

    def _prepare_input_files(self):
        self.logger.debug("Preparing input files for testing purposes...")
        pass

    """
    This helper method was taken from the run.py file, and just prints the output for the user 
    """
    def _print_output(self, result): 
        echo("Printing output...")
        if isinstance(result, types.GeneratorType):
            for r in result:
                if r is not None:
                    echo(json.dumps(r, indent=4))
                else:
                    echo("Something went wrong", fg="red")
        else:
            echo(result)
    

    """
    This helper method checks that the model ID is correct.
    """
    def _check_model_id(self, data): 
        print("Checking model ID...")
        if data["card"]["Identifier"] != self.model_id:
            raise texc.WrongCardIdentifierError(self.model_id)    


    """
    This helper method checks that the slug field is non-empty.
    """
    def _check_model_slug(self, data): 
        print("Checking model slug...")
        if not data["card"]["Slug"]:
            raise texc.EmptyField("slug")
   
    """
    This helper method checks that the description field is non-empty.
    """
    def _check_model_description(self, data): 
        print("Checking model description...")
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
        print("Checking model task...")
        valid_tasks = [[ 'Classification'], [ 'Regression' ], [ 'Generative' ], [ 'Representation' ], 
                       [ 'Similarity' ], [ 'Clustering' ], [ 'Dimensionality reduction' ]]
        if data["card"]["Task"] not in valid_tasks:
            raise texc.InvalidEntry("Task")
    
    """
    This helper method checks that the input field is one of the following valid entries:
        - Compound
        - Protein
        - Text
    """
    def _check_model_input(self, data): 
        print("Checking model input...")
        valid_inputs = [[ 'Compound' ], [ 'Protein' ], [ 'Text' ]]
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
        print("Checking model input shape...")
        valid_input_shapes = ["Single", "Pair", "List", "Pair of Lists", "List of Lists"]
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
        print("Checking model output...")
        valid_outputs = [[ 'Boolean' ], [ 'Compound' ], [ 'Descriptor' ], [ 'Distance' ], [ 'Experimental value' ], 
                         [ 'Image' ], [ 'Other value' ], [ 'Probability' ], [ 'Protein' ], [ 'Score' ], [ 'Text' ]]
        if data["card"]["Output"] not in valid_outputs:
            raise texc.InvalidEntry("Output")

    """
    This helper method checks that the output type is one of the following valid entries:
        - String
        - Float
        - Integer
    """
    def _check_model_output_type(self, data): 
        print("Checking model output type...")
        valid_output_types = [[ 'String' ], [ 'Float' ], [ 'Integer' ]]
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
        print("Checking model output shape...")
        valid_output_shapes = ["Single", "List", "Flexible List", "Matrix", "Serializable Object"]
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
    def check_information(self):
        self.logger.debug("Checking that model information is correct")
        print(BOLD + "Beginning checks for {0} model information:".format(self.model_id) + RESET)
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
        print("SUCCESS! Model information verified.\n")


    """
    Runs the model on a single smiles string and prints the output, or writes it to specified output file
    """
    @throw_ersilia_exception
    def check_single_input(self, output):
        session = Session(config_json=None)
        service_class = session.current_service_class()

        click.echo(BOLD + "Testing model on single smiles input...\n" + RESET)

        input = "COc1ccc2c(NC(=O)Nc3cccc(C(F)(F)F)n3)ccnc2c1"
        mdl = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        result = mdl.run(input=input, output=output, batch_size=100)
        self._print_output(result)

    """
    Generates an example input of 5 smiles using the 'example' command, and then tests the model on that input and prints it
    to the consol.
    """
    @throw_ersilia_exception
    def check_example_input(self, output):
        click.echo(BOLD + "\nTesting model on input of 5 smiles given by 'example' command...\n" + RESET)

        session = Session(config_json=None)
        service_class = session.current_service_class()

        eg = ExampleGenerator(model_id=self.model_id)
        input = eg.example(n_samples=5, file_name=None, simple=True)

        mdl = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        result = mdl.run(input=input, output=output, batch_size=100) 
        self._print_output(result)


    """
    Gets an example input of 5 smiles using the 'example' command, and then runs this same input on the 
    model twice. Then, it checks if the outputs are consistent or not and specifies that to the user. 
    Lastly, it makes sure that the number of outputs equals the number of inputs.  
    """
    @throw_ersilia_exception
    def check_consistent_output(self, output):
        # self.logger.debug("Confirming model produces consistent output...")
        click.echo(BOLD + "\nConfirming model produces consistent output..." + RESET)

        session = Session(config_json=None)
        service_class = session.current_service_class()

        eg = ExampleGenerator(model_id=self.model_id)
        input = eg.example(n_samples=NUM_SAMPLES, file_name=None, simple=True)

        mdl = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        result = mdl.run(input=input, output=output, batch_size=100)
        result2 = mdl.run(input=input, output=output, batch_size=100)

        consistent = True
        zipped = list(zip(result, result2))


        """
        ***
        We can only check for this 5% threshold if the outputs are numbers. If strings, it won't work 
        ***
        """


        for item1, item2 in zipped:
            output1 = item1['output']['mw']
            output2 = item2['output']['mw']
            # check to see if the first and second outputs are within 5% from each other 
            if (100 * (abs(output2 - output1) / ((output1 + output2) / 2)) > DIFFERENCE_THRESHOLD):
                for item1, item2 in zipped:
                    print(item1)
                    print(item2)
                    print('\n')
                raise texc.InconsistentOutputs(self.model_id)

        print("Model output is consistent!")
    
        click.echo(BOLD + "\nConfirming there are same number of outputs as inputs..." + RESET)
        print("Number of inputs:", NUM_SAMPLES)
        print("Number of outputs:", len(zipped))

        if NUM_SAMPLES != len(zipped): 
            raise texc.MissingOutputs()
        else: 
            echo("Number of outputs and inputs are equal!")

    
    @throw_ersilia_exception
    def run_bash(self, output): 
        print("Running the model bash script...")
        # Save current directory - atm, this must be run from root directory (~)
        current_dir = os.getcwd()

        # Create temp directory and clone model 
        with tempfile.TemporaryDirectory() as temp_dir:
            repo_url = 'https://github.com/ersilia-os/{0}.git'.format(self.model_id)  # Replace with the actual GitHub repository URL
            try:
                subprocess.run(['git', 'clone', repo_url, temp_dir], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error while cloning the repository: {e}")

            # Navigate into the temporary directory
            subdirectory_path = os.path.join(temp_dir, "model/framework")
            os.chdir(subdirectory_path)

            # Create temp file
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
                temp_file_path = temp_file.name

            # Run bash script with specified args
            output_path = temp_file_path
            run_path = os.path.join(temp_dir, "model/framework/run.sh")      # path to run.sh
            arg1 = os.path.join(current_dir, "ersilia/test/inputs/compound_singles.csv")      # input
            arg2 = output_path      # output

            try:
                subprocess.run(['bash', run_path, ".", arg1, arg2,], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error while running the bash script: {e}")

            with open(output_path, 'r') as temp_file:
                output_contents = temp_file.read()

            print("Output contents:")
            print(output_contents)

            # maybe instead of printing the contents of the bash, we can run just compare it with an ersilia run


    @throw_ersilia_exception
    def run_using_bash(self): 
        tmp_folder = tempfile.mkdtemp(prefix="eos-")
        run_file = os.path.join(tmp_folder, self.RUN_FILE)
        data_file = os.path.join(tmp_folder, DATA_FILE)
        output_file = os.path.join(tmp_folder, OUTPUT_FILE)

        cur_path = os.path.dirname(os.path.realpath(__file__))
        framework_dir = os.path.join(cur_path, "..", "..", "..", self.model_id, "model", FRAMEWORK_BASEDIR, "run.sh")

        subprocess.run(["chmod +x " + framework_dir, framework_dir, data_file, output_file], shell=True)

    def run(self, output):
        self.check_information()
        self.check_single_input(output)
        self.check_example_input(output)
        self.check_consistent_output(output)
        # self.run_using_bash()
        # self.run_bash(output)


# To do: 
# 1. When it currently prints to an output file, it writes the single output result, then deletes that, then prints the result for the example input. Fix this 
    # to do this, write each output to a temporary file, and then append all of them to a final file and delete the temporary files at the end 
# 2. test it with normal run and then try the bash run.sh, comparing the two outputs 
# 3. Make sure the output matches the expectation of the otuput (ex: if it expects a float, it actually gets a float)
