import csv
import os
import warnings
import tempfile
import traceback
from pathlib import Path
from typing import Callable
import warnings

warnings.filterwarnings("ignore", message="Using slow pure-python SequenceMatcher")

# ruff: noqa
MISSING_PACKAGES = False
try:
    from fuzzywuzzy import fuzz
except ImportError:
    MISSING_PACKAGES = True
# ruff: enable

from .... import throw_ersilia_exception
from ....default import (
    RUN_FILE,
)
from .constants import STATUS_CONFIGS, TABLE_CONFIGS, Checks, TableType
from .setup import SetupService
from .inspect import InspectService
from .io import IOService, PackageInstaller
from .checks import CheckService
from ....default import PREDEFINED_COLUMN_FILE
from ....io.input import ExampleGenerator
from ....utils.exceptions_utils import test_exceptions as texc
from ....utils.spinner import show_loader
from ....utils.terminal import run_command_check_output, run_command
from ....cli import echo

warnings.filterwarnings("ignore", message="Using slow pure-python SequenceMatcher.*")


class RunnerService:
    """
    Service for running model tests and checks.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    logger : logging.Logger
        Logger for logging messages.
    ios_service : IOService
        Instance of IOService for handling input/output operations.
    checkup_service : CheckService
        Instance of CheckService for performing various checks on the model.
    setup_service : SetupService
        Instance of SetupService for setting up the environment and fetching the model repository.
    dir : str
        Directory where the model repository is located.
    from_dir : str
        Directory where the model repository is located.
    model_path : Callable
        Callable to get the model path.
    from_github : bool
        Flag indicating whether to fetch the repository from GitHub.
    from_s3 : bool
        Flag indicating whether to fetch the repository from S3.
    from_dockerhub : bool
        Flag indicating whether to fetch the repository from DockerHub.
    version : str
        Version of the model.
    shallow : bool
        Flag indicating whether to perform shallow checks.
    deep : bool
        Flag indicating whether to perform deep checks.
    report_path : bool
        Flag to specify the path for output as json.
    inspector : InspectService
        Instance of InspectService for inspecting models and their configurations.
    surface : bool
        Flag indicating whether to perform surface checks.
    """

    def __init__(
        self,
        model_id: str,
        logger,
        ios_service: IOService,
        checkup_service: CheckService,
        setup_service: SetupService,
        dir: str,
        from_dir: str,
        model_path: Callable,
        from_github: bool,
        from_s3: bool,
        from_dockerhub: bool,
        version: str,
        shallow: bool,
        deep: bool,
        report_path: str,
        inspector: InspectService,
        surface: bool,
    ):
        self.model_id = model_id
        self.logger = logger
        self.setup_service = setup_service
        self.ios_service = ios_service
        self.console = ios_service.console
        self.checkup_service = checkup_service
        self.model_path = model_path(self.model_id)
        self.dir = dir
        self.from_dir = from_dir
        self.from_github = from_github
        self.from_s3 = from_s3
        self.from_dockerhub = from_dockerhub
        self.version = version
        self.shallow = shallow
        self.deep = deep
        self.report_path = Path(report_path) if report_path else Path().cwd()
        self.report_file = self.report_path / f"{self.model_id}-test.json"
        self.inspector = inspector
        self.example = ExampleGenerator(model_id=self.model_id)
        self.run_using_bash = False
        self.surface = surface
        self.installer = PackageInstaller(self.dir, self.model_id)

    @throw_ersilia_exception()
    def run_model(self, inputs: list, output: str, batch: int):
        """
        Run the model with the given input and output parameters.

        Parameters
        ----------
        input : list
            List of input samples.
        output : str
            Path to the output file.
        batch : int
            Batch size for running the model.

        Returns
        -------
        str
            The output of the command.
        """

        cmd = f"ersilia serve {self.model_id} --no-cache && ersilia run -i '{inputs}' -o {output} -b {str(batch)}"
        out = run_command(cmd)
        return out
    
    @throw_ersilia_exception()
    def delete(self):
        """
        Delete model if existed
        """
        if os.path.exists(self.model_path):
            (
                run_command(
                    [
                        "ersilia",
                        "-v",
                        "delete",
                        self.model_id,
                    ]
                ),
            )

    @throw_ersilia_exception()
    def fetch(self):
        """
        Fetch the model repository from the specified directory.
        """

        def _fetch(model_id, logger):
            loc = (
                ["--from_dir", self.dir]
                if self.from_github or self.from_s3 or self.from_dir
                else ["--from_dockerhub"]
                + (["--version", self.version] if self.version else [])
            )
            self.logger.info(f"Fetching the model from: {loc}")
            run_command(["ersilia", "-v", "fetch", model_id, *loc], quiet=True)

        self.delete()
        _fetch(self.model_id, self.logger)

    def run_example(
        self,
        n_samples: int,
        file_name: str = None,
        simple: bool = True,
        try_predefined: bool = False,
    ):
        """
        Generate example input samples for the model.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate.
        file_name : str, optional
            Name of the file to save the samples.
        simple : bool, optional
            Flag indicating whether to generate simple samples.
        try_predefined : bool, optional
            Flag indicating whether to try predefined samples.

        Returns
        -------
        list
            List of generated input samples.
        """
        examples = self.example.example(
            n_samples=n_samples,
            file_name=file_name,
            simple=simple,
            try_predefined=try_predefined,
        )
        return examples

    @throw_ersilia_exception()
    def run_bash(self):
        """
        Run the model using a bash script and compare the outputs for consistency.

        Raises
        ------
        RuntimeError
            If there is an error during the subprocess execution or output comparison.
        """

        def compute_rmse(y_true, y_pred):
            return sum((yt - yp) ** 2 for yt, yp in zip(y_true, y_pred)) ** 0.5 / len(
                y_true
            )

        def compare_outputs(bsh_data, ers_data):
            _completed_status, _rmse = [], []
            columns = set(bsh_data[0].keys()) & set(ers_data[0].keys())
            for column in columns:
                bv = [row[column] for row in bsh_data]
                ev = [row[column] for row in ers_data]

                if all(isinstance(val, (int, float)) for val in bv + ev):
                    rmse = compute_rmse(bv, ev)
                    _rmse.append(rmse)

                    if rmse > 0.1:
                        rmse_perc = round(rmse * 100, 2)
                        _completed_status.append(
                            (
                                f"RMSE-{column}",
                                f"RMSE > 10%{rmse_perc}%",
                                str(STATUS_CONFIGS.FAILED),
                            )
                        )
                        raise texc.InconsistentOutputs(self.model_id)

                elif all(isinstance(val, str) for val in bv + ev):
                    if not all(
                        self._compare_string_similarity(a, b, 95)
                        for a, b in zip(bv, ev)
                    ):
                        _completed_status.append(
                            ("String Similarity", "< 95%", str(STATUS_CONFIGS.FAILED))
                        )
                        raise texc.InconsistentOutputs(self.model_id)
                    _completed_status.append(
                        ("String Similarity", "> 95%", str(STATUS_CONFIGS.PASSED))
                    )

            rmse = sum(_rmse) / len(_rmse) if _rmse else 0
            rmse_perc = round(rmse * 100, 2)
            _completed_status.append(
                ("RMSE-MEAN", f"RMSE < 10% | {rmse_perc}%", str(STATUS_CONFIGS.PASSED))
            )

            return _completed_status

        def read_csv(path, flag=False):
            if not os.path.exists(path):
                raise FileNotFoundError(
                    f"File not found this might be due to error happened when executing the run.sh: {path}"
                )
            try:
                with open(path, "r") as file:
                    self.logger.info("Reading the lines")
                    lines = file.readlines()

                if not lines:
                    self.logger.error(f"File at {path} is empty.")
                    return [], "File is empty"

                headers = lines[0].strip().split(",")
                if flag:
                    headers = headers[2:]

                data = []

                for line in lines[1:]:
                    values = line.strip().split(",")
                    values = values[2:] if flag else values

                    def parse(value):
                        try:
                            v = int(value)
                            return v
                        except ValueError:
                            pass

                        try:
                            v = float(value)
                            return v
                        except ValueError:
                            pass

                        if isinstance(value, str):
                            return value

                    try:
                        _values = [parse(x) for x in values]
                    except ValueError as e:
                        return [], f"Invalid value detected in CSV file: {e}"

                    data.append(dict(zip(headers, _values)))

                return data, None

            except Exception as e:
                raise RuntimeError(f"Failed to read CSV from {path}.") from e

        def read_logs(path):
            if not os.path.exists(path):
                raise FileNotFoundError(f"File not found: {path}")
            with open(path, "r") as file:
                return file.readlines()

        def rename_col(path, old_col_name="input", new_col_name="smiles"):
            with open(path, "r", newline="") as infile:
                reader = csv.reader(infile)
                rows = list(reader)
                self.logger.info(f"Rows before rename: {rows}")

                if rows and old_col_name in rows[0]:
                    col_index = rows[0].index(old_col_name)
                    rows[0][col_index] = new_col_name

                with open(path, "w", newline="") as outfile:
                    writer = csv.writer(outfile)
                    writer.writerows(rows)

        with tempfile.TemporaryDirectory() as temp_dir:
            model_path = os.path.join(self.dir)
            output_path = os.path.join(temp_dir, "ersilia_output.csv")
            output_log_path = os.path.join(temp_dir, "output.txt")
            error_log_path = os.path.join(temp_dir, "error.txt")
            input_file_path = os.path.join(temp_dir, "example_file.csv")
            temp_script_path = os.path.join(temp_dir, "script.sh")
            bash_output_path = os.path.join(temp_dir, "bash_output.csv")

            input_file_path = IOService._get_input_file_path(self.dir)
            rename_col(input_file_path)

            run_sh_path = os.path.join(model_path, "model", "framework", RUN_FILE)
            if not os.path.exists(run_sh_path):
                self.logger.warning(
                    f"{RUN_FILE} not found at {run_sh_path}. Skipping bash run."
                )
                return
            
            self.installer._install_packages_from_dir()

            bash_script = f"""
                source {self._conda_prefix(self._is_base())}/etc/profile.d/conda.sh
                conda activate {self.model_id}
                cd {os.path.dirname(run_sh_path)}
                bash run.sh . {input_file_path} {bash_output_path} > {output_log_path} 2> {error_log_path}
                conda deactivate
                """

            with open(temp_script_path, "w") as script_file:
                script_file.write(bash_script)

            self.logger.debug(f"\nRunning bash script: {temp_script_path}\n")
            out = run_command(["bash", temp_script_path])
            self.logger.info(f"Bash script subprocess output: {out}")
            logs = read_logs(error_log_path)
            formatted_error = "".join(logs)
            if formatted_error:
                echo(
                    f"Running bash: {formatted_error}",
                    fg="yellow",
                    bold=True,
                )
            bsh_data, _ = read_csv(bash_output_path)
            self.logger.debug("Running model for bash data consistency checking")
            cmd = f"ersilia serve {self.model_id} --no-cache && ersilia -v run -i '{input_file_path}' -o {output_path}"
            out = run_command(cmd)
            ers_data, _ = read_csv(output_path, flag=True)
            self.checkup_service.original_smiles_list = (
                self.checkup_service._get_original_smiles_list("csv", input_file_path)
            )
            check_status = self.checkup_service._check_csv(
                output_path, input_type="csv"
            )
            if check_status[-1] == str(STATUS_CONFIGS.FAILED):
                return [
                    (
                        (
                            Checks.RUN_BASH.value,
                            check_status[1],
                            str(STATUS_CONFIGS.FAILED),
                        )
                    )
                ]
            status = compare_outputs(bsh_data, ers_data)
            return status

    @staticmethod
    def _default_env():
        if "CONDA_DEFAULT_ENV" in os.environ:
            return os.environ["CONDA_DEFAULT_ENV"]
        return None

    @staticmethod
    def _conda_prefix(is_base):
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

    def _is_base(self):
        default_env = self._default_env()
        self.logger.debug(f"Default environment: {default_env}")
        return default_env == "base"

    def _compare_string_similarity(self, str1, str2, threshold):
        similarity = fuzz.ratio(str1, str2)
        return similarity >= threshold

    def run(self):
        """
        Run the model tests and checks.
        """
        results = []
        try:
            self._configure_environment()
            self.setup_service.get_model()

            basic_results = self._perform_basic_checks()
            results.extend(basic_results)

            if self.surface:
                surface_results = self._perform_surface_check()
                results.extend(surface_results)

            if self.shallow:
                results.extend(self._perform_surface_check())
                results.extend(self._perform_shallow_checks())

            if self.deep:
                results.extend(self._perform_surface_check())
                results.extend(self._perform_shallow_checks())
                deep_result = self._perform_deep_checks()
                results.append(deep_result)

            self.ios_service.collect_and_save_json(results, self.report_file)
            echo("Model tests and checks completed.", fg="green", bold=True)
            
            echo("Deleting model...", fg="yellow", bold=True)
            self.delete()
            echo("Model successfully deleted", fg="green", bold=True)

        except Exception as error:
            tb = traceback.format_exc()
            error_info = {"exception": str(error), "traceback": tb}
            results.append(error_info)
            echo(f"An error occurred: {error}\nTraceback:\n{tb}", fg="red", bold=True)
            echo("Deleting model...", fg="yellow", bold=True)
            self.delete()
            echo("Model successfully deleted", fg="green", bold=True)

            self.ios_service.collect_and_save_json(results, self.report_file)

    def _configure_environment(self):
        if self.from_dockerhub:
            self.setup_service.from_github = True

    def _perform_basic_checks(self):
        results = []

        self.checkup_service.check_information()
        results.append(
            self._generate_table_from_check(
                TableType.MODEL_INFORMATION_CHECKS, self.ios_service.check_results
            )
        )
        self.ios_service.check_results.clear()

        self.checkup_service.check_files()
        results.append(
            self._generate_table_from_check(
                TableType.MODEL_FILE_CHECKS, self.ios_service.check_results
            )
        )

        results.append(self._docker_yml_column_name_check())
        results.append(self._log_directory_sizes())
        return results

    @show_loader(text="Performing surface checks", color="cyan")
    def _perform_surface_check(self):
        self.fetch()
        results = []

        if self.from_github or self.from_s3 or self.from_dir:
            env_result = self._log_env_sizes()
            results.append(
                self._generate_table_from_check(TableType.MODEL_SIZES, env_result)
            )

        if self.from_dockerhub:
            tag = self.version if self.version else "latest"
            docker_size, message = self.ios_service.calculate_image_size(tag=tag)
            results.append(
                self._generate_table_from_check(
                    TableType.MODEL_SIZES,
                    [(Checks.IMAGE_SIZE.value, message, docker_size)],
                )
            )

        simple_output = self.checkup_service.check_simple_model_output(self.run_model)
        results.append(
            self._generate_table_from_check(TableType.MODEL_RUN_CHECK, simple_output)
        )

        return results

    @show_loader(text="Performing shallow checks", color="cyan")
    def _perform_shallow_checks(self):
        results = []

        model_output = self.checkup_service.check_model_output_content(
            self.run_example, self.run_model
        )
        results.append(
            self._generate_table_from_check(TableType.MODEL_OUTPUT, model_output)
        )
        validations = []
        if "Fixed" in self.ios_service.get_output_consistency():
            res = self._run_single_and_example_input_checks()

            validations.append(
                self._generate_table_from_check(TableType.SHALLOW_CHECK_SUMMARY, res)
            )

        bash_results = self.run_bash()
        validations.append(
            self._generate_table_from_check(TableType.CONSISTENCY_BASH, bash_results)
        )
        results.extend(validations)
        return results

    @show_loader(text="Performing deep checks", color="cyan")
    def _perform_deep_checks(self):
        performance_data = self.inspector.run(["computational_performance_tracking"])
        return self._generate_table_from_check(
            TableType.COMPUTATIONAL_PERFORMANCE, performance_data, large=True
        )

    def _docker_yml_column_name_check(self):
        column_file_path = os.path.join(self.dir, PREDEFINED_COLUMN_FILE)
        example_output_file_path = IOService._get_output_file_path(self.dir)
        data = self.checkup_service.compare_csv_columns(
            column_csv=column_file_path, csv_file=example_output_file_path
        )
        docker_check_data = self.inspector.run(["docker_check"]).pop(0)
        docker_check_data = [
            (docker_check_data[0], Checks.DEPENDENCY_PINNED.value, docker_check_data[1])
        ]
        data.extend(docker_check_data)
        return self._generate_table_from_check(TableType.DEPENDECY_COLUMN_CHECK, data)

    def _log_env_sizes(self):
        env_size = self.ios_service.get_env_sizes()
        return [(Checks.ENV_SIZE.value, Checks.SIZE_CACL_SUCCESS.value, env_size)]

    def _log_directory_sizes(self):
        directory_size = self.ios_service.get_directories_sizes()
        return self._generate_table_from_check(
            TableType.MODEL_DIRECTORY_SIZES,
            [(Checks.DIR_SIZE.value, Checks.TOTAL_DIR_SIZE.value, directory_size)],
        )

    def _run_single_and_example_input_checks(self):
        results = []
        results.extend(
            self.checkup_service.check_consistent_output(
                self.run_example, self.run_model
            )
        )
        return results

    def _generate_table_from_check(self, table_type, rows, large=False):
        config = TABLE_CONFIGS[table_type].__dict__
        return self.ios_service._generate_table(
            title=config["title"],
            headers=config["headers"],
            rows=rows,
            large_table=large,
        )
