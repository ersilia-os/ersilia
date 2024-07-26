import os

from ....utils.terminal import run_command_check_output, run_command
from ....utils.conda import SimpleConda
from ....utils.exceptions_utils.fetch_exceptions import StandardModelExampleError

from ....default import (
    EXAMPLE_STANDARD_INPUT_CSV_FILENAME,
    EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME,
)

from .... import ErsiliaBase
from .... import throw_ersilia_exception


class ModelStandardExample(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    @throw_ersilia_exception
    def _check_file_exists(self, output_csv):
        if not os.path.exists(output_csv):
            raise StandardModelExampleError(
                model_id=self.model_id, file_name=output_csv
            )
        self.logger.debug("File {0} created successfully!".format(output_csv))

    def run(self):
        self.logger.debug("Running standard CSV example")
        path = self._model_path(model_id=self.model_id)
        input_csv = os.path.join(path, EXAMPLE_STANDARD_INPUT_CSV_FILENAME)
        output_csv = os.path.join(path, EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME)
        run_log = os.path.join(path, "standard_run.log")
        self.logger.debug(input_csv)
        self.logger.debug(output_csv)
        commands = [
            "ersilia serve {0}".format(self.model_id),
            "ersilia example -n 3 -c -f {0}".format(input_csv),
            "ersilia -v run -i {0} -o {1} > {2} 2>&1".format(
                input_csv, output_csv, run_log
            ),
            "ersilia close",
        ]
        cmd_output = run_command_check_output("ersilia --help")
        self.logger.debug(cmd_output)
        if "Welcome to Ersilia" in cmd_output:
            self.logger.debug("No need to use Conda!")
            cmd = " && ".join(commands)
            run_command(cmd)
        else:
            self.logger.debug("Will run this through Conda")
            env_name = os.environ.get("CONDA_DEFAULT_ENV")
            self.logger.debug("The environment name is {0}".format(env_name))
            SimpleConda().run_commandlines(env_name, commands)
        self._check_file_exists(output_csv=output_csv)
        self.logger.debug("Removing log file: {0}".format(run_log))
        os.remove(run_log)
