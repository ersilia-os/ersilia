import os
import json
import click
import tempfile
import types
from ..cli import echo

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
    def check_single_input(self):
        self.logger.debug("Testing model on the following single smiles input:  COc1ccc2c(NC(=O)Nc3cccc(C(F)(F)F)n3)ccnc2c1")
        click.echo("Testing model on the following single smiles input:  COc1ccc2c(NC(=O)Nc3cccc(C(F)(F)F)n3)ccnc2c1...")

        session = Session(config_json=None)
        service_class = session.current_service_class()

        input = "COc1ccc2c(NC(=O)Nc3cccc(C(F)(F)F)n3)ccnc2c1"
        mdl = ErsiliaModel(self.model_id, service_class=service_class, config_json=None)
        result = mdl.run(input=input, output=None, batch_size=100)

        if isinstance(result, types.GeneratorType):
            for result in mdl.run(input=input, output=None, batch_size=100):
                if result is not None:
                    echo(json.dumps(result, indent=4))
                else:
                    echo("Something went wrong", fg="red")
        else:
            echo(result)


    def run(self):
        self.check_information()
        self.check_single_input()
        print('all tested!')
