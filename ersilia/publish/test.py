import os
import json
import click
import tempfile
import types
from ..cli import echo

from ..io.input import ExampleGenerator
from .. import ErsiliaBase
from .. import throw_ersilia_exception
from .. import ErsiliaModel

from ..utils.exceptions_utils import test_exceptions as texc
from ..core.session import Session

from ..default import INFORMATION_FILE


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
        with open(json_file, "r") as f:
            data = json.load(f)
        return data

    def _prepare_input_files(self):
        self.logger.debug("Preparing input files for testing purposes...")
        pass

    @throw_ersilia_exception
    def check_information(self):
        self.logger.debug("Checking that model information is correct")
        json_file = os.path.join(self._dest_dir, self.model_id, INFORMATION_FILE)
        with open(json_file, "r") as f:
            data = json.load(f)
        if data["card"]["Identifier"] != self.model_id:
            raise texc.WrongCardIdentifierError(self.model_id)
        return data

    @throw_ersilia_exception
    def check_single_input(self, output):
        # self.logger.debug("Testing model on custom example input with 10 smiles...")
        click.echo("Testing model on custom example input with 5 smiles...")

        session = Session(config_json=None)
        service_class = session.current_service_class()

        eg = ExampleGenerator(model_id=self.model_id)
        input = eg.example(n_samples=5, file_name=None, simple=True)

        mdl = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        result = mdl.run(input=input, output=output, batch_size=100)

        if isinstance(result, types.GeneratorType):
            for result in mdl.run(input=input, output=output, batch_size=100):
                if result is not None:
                    echo(json.dumps(result, indent=4))
                else:
                    echo("Something went wrong", fg="red")
            # print('all tested!')
        else:
            echo(result)

            # need to adjust it to actually put the output into the output file if specified

            print('all tested!')


    def run(self, output):
        self.check_information()
        self.check_single_input(output)
