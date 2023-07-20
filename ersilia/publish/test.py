import os
import json
import click
import tempfile
import types
import time
from ..cli import echo

from ..io.input import ExampleGenerator
from .. import ErsiliaBase
from .. import throw_ersilia_exception
from .. import ErsiliaModel

from ..utils.exceptions_utils import test_exceptions as texc
from ..core.session import Session

from ..default import INFORMATION_FILE
NUM_SAMPLES = 5

class ModelTester(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        self._info = self._read_information()
        self._input = self._info["card"]["Input"]
        self._prepare_input_files()

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
        if isinstance(result, types.GeneratorType):
            for r in result:
                if r is not None:
                    echo(json.dumps(r, indent=4))
                else:
                    echo("Something went wrong", fg="red")
        else:
            echo(result)


    """
    Check the model information to make sure it's correct. Need to add more here.
    """
    @throw_ersilia_exception
    def check_information(self):
        self.logger.debug("Checking that model information is correct")
        json_file = os.path.join(self._dest_dir, self.model_id, INFORMATION_FILE)
        with open(json_file, "r") as f:
            data = json.load(f)
        if data["card"]["Identifier"] != self.model_id:
            raise texc.WrongCardIdentifierError(self.model_id)
        return data

    """
    Runs the model on a single smiles string and prints the output, or writes it to specified output file
    """
    @throw_ersilia_exception
    def check_single_input(self, output):
        session = Session(config_json=None)
        service_class = session.current_service_class()

        click.echo("Testing model on single smiles input...\n")

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
        click.echo("\nTesting model on input of 5 smiles given by 'example' command...\n")

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
        click.echo("\nConfirming model produces consistent output...")

        session = Session(config_json=None)
        service_class = session.current_service_class()

        eg = ExampleGenerator(model_id=self.model_id)
        input = eg.example(n_samples=NUM_SAMPLES, file_name=None, simple=True)

        mdl = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        result = mdl.run(input=input, output=output, batch_size=100)
        result2 = mdl.run(input=input, output=output, batch_size=100)

        consistent = True
        zipped = list(zip(result, result2))

        for item1, item2 in zipped:
            if (item1 != item2):
                consistent = False
                break

        if consistent:
            print("Model output is consistent!")
        else: 
            print("Model output is inconsistent. Please review the output shown below to see if there is an issue with the model. If the differences are small, there may be no issue!\n")
            for item1, item2 in zipped:
                print(item1)
                print(item2)
                print('\n')

        click.echo("\nConfirming there are same number of outputs as inputs...")
        print("Number of inputs:", NUM_SAMPLES)
        print("Number of outputs:", len(zipped))

        if NUM_SAMPLES != len(zipped): 
            raise texc.MissingOutputs(self.model_id)
        else: 
            echo("Number of outputs and inputs are equal!")


    @throw_ersilia_exception
    def run_bash(self, output): 
        pass


    def run(self, output):
        self.check_information()
        self.check_single_input(output)
        self.check_example_input(output)
        self.check_consistent_output(output)
        # self.run_bash(output)


# To do: 
# 1. When it currently prints to an output file, it writes the single output result, then deletes that, then prints the result for the example input. Fix this
# 2. test it with normal run and then try the bash run.sh, comparing the two outputs 
# 3. speed it up by making it so the check_consistent_output function doesn't have to re-generate the example input and then run the model on that input twice
# 4. for each major aspect of the test module, I'd like the header when showing the user what it's doing to be more like a header (bold, different color, etc.)
