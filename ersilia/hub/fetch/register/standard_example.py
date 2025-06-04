import csv
import os

from .... import ErsiliaBase, throw_ersilia_exception
from ....default import (
    EXAMPLE_STANDARD_INPUT_CSV_FILENAME,
    EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME,
)
from ....utils.conda import SimpleConda
from ....utils.exceptions_utils.fetch_exceptions import StandardModelExampleError
from ....utils.terminal import run_command, run_command_check_output


class ModelStandardExample(ErsiliaBase):
    """
    ModelStandardExample is responsible for running standard CSV examples for models.

    Parameters
    ----------
    model_id : str
        The ID of the model to run the example for.
    config_json : dict
        Configuration settings for the example runner.

    Examples
    --------
    .. code-block:: python

        example_runner = ModelStandardExample(
            model_id="model123", config_json=config
        )
        example_runner.run()
    """

    def __init__(self, model_id: str, config_json: dict):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def _less_on_file(self, output_csv: str) -> list:
        max_rows = 5
        max_cols = 10
        lines = []
        with open(output_csv, newline="", encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile)
            data = []
            for i, r in enumerate(reader):
                data += [r]
                if i > max_rows:
                    break
            col_widths = [max(len(str(item)) for item in col) for col in zip(*data)]
            for row in data:
                if len(row[:max_cols]) > len(col_widths):
                    max_cols = len(col_widths)
                lines += [
                    " | ".join(
                        f"{int(item):<{col_widths[i]}}"
                        if item.isdigit()
                        else f"{float(item):<{col_widths[i]}.2f}"
                        if item.replace(".", "", 1).isdigit()
                        else f"{item:<{col_widths[i]}}"
                        for i, item in enumerate(row[:max_cols])
                    )
                ]
        return lines

    @throw_ersilia_exception(exit=False)
    def _check_file_exists(self, output_csv: str):
        if not os.path.exists(output_csv):
            raise StandardModelExampleError(
                model_id=self.model_id, file_name=output_csv
            )
        self.logger.debug("File {0} created successfully!".format(output_csv))
        self.logger.debug("This is the output (maximum 5 rows, 10 columns)")
        lines = self._less_on_file(output_csv=output_csv)
        for l in lines:
            self.logger.debug(l)

    def run(self):
        """
        Run the standard CSV example for the model.

        This method runs the standard CSV example for the model, generating input and output CSV files.
        """
        self.logger.debug("Running standard CSV example")
        path = self._model_path(model_id=self.model_id)
        input_csv = os.path.join(path, EXAMPLE_STANDARD_INPUT_CSV_FILENAME)
        output_csv = os.path.join(path, EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME)
        run_log = os.path.join(path, "standard_run.log")
        self.logger.debug(input_csv)
        self.logger.debug(output_csv)
        commands = [
            "ersilia serve {0} --disable-local-cache".format(self.model_id),
            "ersilia example -n 3 -f {0} --simple".format(input_csv),
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

        self.logger.info(f"Run log: {open(run_log).read()}")
        self._check_file_exists(output_csv=output_csv)
        self.logger.debug("Removing log file: {0}".format(run_log))
        os.remove(run_log)
