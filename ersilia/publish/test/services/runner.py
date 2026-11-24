import csv
import filecmp
import os
import sys
import warnings
import tempfile
import traceback
from pathlib import Path
from typing import Callable

warnings.filterwarnings("ignore", message="Using slow pure-python SequenceMatcher")

# ruff: noqa
MISSING_PACKAGES = False
try:
    from fuzzywuzzy import fuzz
except ImportError:
    MISSING_PACKAGES = True
# ruff: enable

from ....default import (
    RUN_FILE,
)
from .constants import STATUS_CONFIGS, TABLE_CONFIGS, Checks, TableType
from .setup import SetupService
from .inspect import InspectService
from .io import IOService, PackageInstaller
from .checks import CheckService
from ....default import (
    EXAMPLE_STANDARD_INPUT_CSV_FILENAME,
    EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME,
    PREDEFINED_COLUMN_FILE,
    _CONDA_BOOTSTRAP,
)
from ....io.input import ExampleGenerator
from ....utils.terminal import run_command_check_output, run_command
from ....cli import echo
from ....store.utils import echo_exceptions, ClickInterface

warnings.filterwarnings("ignore", message="Using slow pure-python SequenceMatcher.*")

_NON_FINITE_STRINGS = ("nan", "+nan", "-nan", "none", "inf", "+inf", "-inf")


def _fetched_example_output(model_path, test_input):
    """
    Find the output that 'ersilia fetch' produced for the example input.

    Fetch runs the model on its copy of ``run_input.csv`` and keeps the
    result. The test can reuse it instead of running the model again, but
    only when fetch ran on exactly the input the test would use: the copy
    inside a Docker image can differ from the repository's.

    Parameters
    ----------
    model_path : str
        The fetched model's folder (``~/eos/dest/<model_id>``).
    test_input : str or None
        The example input the test would run (the repository's
        ``run_input.csv``).

    Returns
    -------
    str or None
        The path to fetch's output, or None if it cannot be reused.
    """
    fetch_input = os.path.join(model_path, EXAMPLE_STANDARD_INPUT_CSV_FILENAME)
    fetch_output = os.path.join(model_path, EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME)
    if not test_input or not os.path.isfile(test_input):
        return None
    if not os.path.isfile(fetch_input) or not os.path.isfile(fetch_output):
        return None
    if os.path.getsize(fetch_output) == 0:
        return None
    if not filecmp.cmp(fetch_input, test_input, shallow=False):
        return None
    return fetch_output


def _parse_cell(value: str):
    """
    Convert a CSV cell to int or float when it looks numeric, otherwise keep the string.

    Empty cells and non-finite markers such as ``nan`` or ``inf`` stay as strings.

    Parameters
    ----------
    value : str
        Raw cell value.

    Returns
    -------
    int, float or str
        The parsed value.
    """
    if value == "":
        return value
    try:
        i = int(value)
        if str(i) == value or (
            value.startswith(("+", "-")) and str(i) == value.lstrip("+")
        ):
            return i
    except ValueError:
        pass
    try:
        f = float(value)
        if value.lower() not in _NON_FINITE_STRINGS:
            return f
    except ValueError:
        pass
    return value


def _read_output_csv(path: str, skip_columns: int = 0):
    """
    Read a model output CSV into a list of row dictionaries with parsed values.

    Uses the csv module, so a quoted cell containing commas (for example a
    CXSMILES such as ``CCO |&1:8,10|``) stays in a single column.

    Parameters
    ----------
    path : str
        Path to the CSV file.
    skip_columns : int, optional
        Number of leading columns to drop, e.g. 2 for ``key`` and ``input`` in
        Ersilia output (default = 0).

    Returns
    -------
    tuple
        ``(rows, error)``, where ``rows`` is a list of dicts keyed by column name
        and ``error`` is None, or a message if the file is empty.
    """
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header:
            return [], "File is empty"
        header = header[skip_columns:]
        rows = []
        for values in reader:
            if not values:
                continue
            values = values[skip_columns:]
            rows.append(dict(zip(header, (_parse_cell(v) for v in values))))
    return rows, None


def _normalize_quotes(obj):
    """
    Strip one pair of surrounding single or double quotes from strings.

    Parameters
    ----------
    obj : str or list
        A string, or a list of strings and numbers.

    Returns
    -------
    str or list
        The input with quotes stripped from each string; numbers are unchanged.
    """

    def strip_quotes(s):
        if isinstance(s, (int, float)):
            return s
        if (s.startswith('"') and s.endswith('"')) or (
            s.startswith("'") and s.endswith("'")
        ):
            return s[1:-1]
        return s

    if isinstance(obj, list):
        return [strip_quotes(s) for s in obj]
    elif isinstance(obj, str):
        return strip_quotes(obj)
    else:
        raise TypeError("Input must be a string or a list of strings")


def _compute_rmse(y_true, y_pred):
    return sum((yt - yp) ** 2 for yt, yp in zip(y_true, y_pred)) ** 0.5 / len(y_true)


def _compare_outputs(bsh_data, ers_data, is_similar: Callable):
    """
    Compare bash and Ersilia run outputs column by column.

    Numeric columns are compared by RMSE and text columns by string similarity.
    Cells that are empty in both outputs are left out. A column that mixes numbers
    and text fails instead of being skipped, because that usually means the CSV
    was not read correctly.

    Parameters
    ----------
    bsh_data : list of dict
        Rows from the bash run output.
    ers_data : list of dict
        Rows from the Ersilia run output.
    is_similar : Callable
        ``is_similar(a, b)`` returns True if two strings count as the same output.

    Returns
    -------
    list of tuple
        Status rows ``(name, detail, status)`` for the test report.
    """
    failed = str(STATUS_CONFIGS.FAILED)
    if not bsh_data or not ers_data:
        msg = (
            f"Cannot compare outputs: bash output has {len(bsh_data)} rows, "
            f"ersilia output has {len(ers_data)} rows"
        )
        echo_exceptions(msg, ClickInterface())
        return [(Checks.RUN_BASH.value, msg, failed)]

    columns = sorted(set(bsh_data[0].keys()) & set(ers_data[0].keys()))
    if not columns:
        msg = "Cannot compare outputs: bash and ersilia outputs share no column names"
        echo_exceptions(msg, ClickInterface())
        return [(Checks.RUN_BASH.value, msg, failed)]

    _completed_status, _rmse = [], []
    for column in columns:
        bv = [row.get(column, "") for row in bsh_data]
        ev = _normalize_quotes([row.get(column, "") for row in ers_data])

        for idx, (b_val, e_val) in enumerate(zip(bv, ev), start=1):
            if type(b_val) is not type(e_val):
                if type(b_val) is int and type(e_val) is float:
                    continue
                msg = (
                    f"Datatype mismatch for column '{column}' at row {idx}: "
                    f"bash value={b_val!r} ({type(b_val).__name__}), "
                    f"ersilia value={e_val!r} ({type(e_val).__name__})"
                )
                echo_exceptions(msg, ClickInterface())
                return [(f"Type Check-{column}", msg, failed)]

        # Rows where both outputs are empty (e.g. variable-length outputs) carry
        # nothing to compare; a one-sided empty cell already failed the type check
        # above for numeric columns, and is compared as text for string columns.
        pairs = [(b, e) for b, e in zip(bv, ev) if not (b == "" and e == "")]
        values = [v for pair in pairs for v in pair]
        if not values:
            continue

        if all(isinstance(v, (int, float)) for v in values):
            rmse = _compute_rmse([b for b, _ in pairs], [e for _, e in pairs])
            _rmse.append(rmse)

            if rmse > 0.1:
                rmse_perc = round(rmse * 100, 2)
                _completed_status.append(
                    (f"RMSE-{column}", f"RMSE > 10%: {rmse_perc}%", failed)
                )
                echo_exceptions(
                    "Model output is inconsistent between bash and ersilia. Skipped the checks!",
                    ClickInterface(),
                )
                return _completed_status

        elif all(isinstance(v, str) for v in values):
            if not all(is_similar(b, e) for b, e in pairs):
                _completed_status.append(("String Similarity", "< 95%", failed))
                echo_exceptions(
                    "Model output is inconsistent between bash and ersilia. Skipped the checks!",
                    ClickInterface(),
                )
                return _completed_status

            _completed_status.append(
                ("String Similarity", "> 95%", str(STATUS_CONFIGS.PASSED))
            )

        else:
            n_text = sum(isinstance(v, str) for v in values)
            msg = (
                f"Column '{column}' mixes numbers and text "
                f"({len(values) - n_text} numeric, {n_text} text values) and cannot "
                "be compared. This usually means the output CSV was not parsed correctly."
            )
            echo_exceptions(msg, ClickInterface())
            return [(f"Type Check-{column}", msg, failed)]

    mean_rmse = sum(_rmse) / len(_rmse) if _rmse else 0
    mean_rmse_perc = round(mean_rmse * 100, 2)
    _completed_status.append(
        (
            "RMSE-MEAN",
            f"RMSE < 10% | {mean_rmse_perc}%",
            str(STATUS_CONFIGS.PASSED),
        )
    )
    return _completed_status


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
    inspect : bool
        Flag indicating whether to perform inspect checks.
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
        inspect: bool,
        permissive: bool = False,
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
        self.inspect = inspect
        self.shallow = shallow
        self.deep = deep
        self.report_path = Path(report_path) if report_path else Path().cwd()
        self.report_file = self.report_path / f"{self.model_id}-test.json"
        self.inspector = inspector
        self.example = ExampleGenerator(model_id=self.model_id)
        self.run_using_bash = False
        self.fetch_output = None
        self.surface = surface
        self.permissive = permissive
        self.installer = PackageInstaller(self.dir, self.model_id)

    def serve_model(self):
        """
        Serve the model and return the command output.

        Returns
        -------
        str
            The output of the command.
        """

        cmd = f"ersilia serve {self.model_id} --disable-cache"
        out = run_command(cmd)
        return out

    def run_model(self, inputs: str, output: str, batch: int):
        """
        Run the model with the given input and output parameters.

        Parameters
        ----------
        input : str
            A path to input file
        output : str
            Path to the output file.
        batch : int
            Batch size for running the model.

        Returns
        -------
        CommandResult
            Named tuple (returncode, stdout, stderr) of the command. The return
            code is that of ``ersilia run`` when serving succeeded, e.g.
            NO_INPUT_PROCESSED_EXIT_CODE when no input gave a result.
        """

        cmd = f"ersilia serve {self.model_id} --disable-cache && ersilia run -i '{inputs}' -o {output} -b {str(batch)}"
        out = run_command(cmd)
        return out

    def _reusable_fetch_output(self):
        """
        Get the example output from fetch, if the test can reuse it.

        Returns
        -------
        str or None
            The path to the output, or None if the model must run again.
        """
        test_input = IOService._get_input_file_path(self.dir) if self.dir else None
        path = _fetched_example_output(self.model_path, test_input)
        if path is None:
            self.logger.info(
                "The example output from fetch cannot be reused; the model will run again."
            )
        return path

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
            cmd = " ".join(["ersilia", "-v", "fetch", model_id, *loc])
            self.logger.debug(f"Running fetch command for testing: {cmd}")
            out = run_command(cmd, quiet=True)
            return out

        self.delete()
        out = _fetch(self.model_id, self.logger)
        return out

    def run_example(
        self,
        n_samples: int,
        file_name: str = None,
        simple: bool = True,
        mode: str = "random",
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
        mode : str, optional
            Mode for generating samples. Can be "predefined", "deterministic", or "random" (default = "random").

        Returns
        -------
        list
            List of generated input samples.
        """
        examples = self.example.example(
            n_samples=n_samples,
            file_name=file_name,
            mode=mode,
        )
        return examples

    def run_bash(self):
        """
        Run the model using a bash script and compare the outputs for consistency.

        Raises
        ------
        RuntimeError
            If there is an error during the subprocess execution or output comparison.
        """

        def read_csv(path, flag=False):
            if not os.path.exists(path):
                msg = f"File not found this might be due to error happened when executing the run.sh: {path}"
                echo_exceptions(msg, ClickInterface())
                return [], msg
            try:
                self.logger.info(f"Reading output CSV: {path}")
                data, error = _read_output_csv(path, skip_columns=2 if flag else 0)
            except Exception as e:
                msg = f"Failed to read CSV from {path}: {e}"
                echo_exceptions(msg, ClickInterface())
                return [], msg
            if error:
                self.logger.error(f"File at {path}: {error}")
            return data, error

        def read_logs(path):
            if not os.path.exists(path):
                return "No error detected!"
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

            run_sh_path = os.path.abspath(
                os.path.join(model_path, "model", "framework", RUN_FILE)
            )
            input_file_path = os.path.abspath(input_file_path)
            if not os.path.exists(run_sh_path):
                self.logger.warning(
                    f"{RUN_FILE} not found at {run_sh_path}. Skipping bash run."
                )
                return

            self.installer.install_packages_from_dir()

            self.logger.debug("The self.dir is: {0}".format(self.dir))
            self.logger.debug("Input file path: {0}".format(input_file_path))
            self.logger.debug("Run script path: {0}".format(run_sh_path))
            self.logger.debug("Output path: {0}".format(output_path))
            bash_script = f"""#!/usr/bin/env bash
            set -euo pipefail

            log_out="{output_log_path}"
            log_err="{error_log_path}"
            mkdir -p "$(dirname "$log_out")" "$(dirname "$log_err")"
            : >"$log_out"
            : >"$log_err"
            exec > >(tee -a "$log_out") 2> >(tee -a "$log_err" >&2)

            ENV_NAME="{self.model_id}"
            RUN_DIR="{os.path.dirname(run_sh_path)}"
            INPUT_FILE="{input_file_path}"
            OUTPUT_FILE="{bash_output_path}"

            echo "[INFO] Runner arch: $(uname -m)"
            echo "[INFO] Target conda env: $ENV_NAME"
            echo "[INFO] Shell: ${{SHELL:-'(unknown)'}}"
            echo "[INFO] PWD: $(pwd)"
            echo "[INFO] Initial conda on PATH: $(command -v conda || echo '(not found)')"
            echo "[INFO] Initial CONDA_SHLVL: ${{CONDA_SHLVL:-0}}"
            echo "[INFO] Initial CONDA_PREFIX: ${{CONDA_PREFIX:-'(unset)'}}"
            echo "[INFO] Initial CONDA_PREFIX_1: ${{CONDA_PREFIX_1:-'(unset)'}}"

            {_CONDA_BOOTSTRAP}

            export CONDA_NO_PLUGINS=true

            echo "[INFO] Conda resolved to: $(command -v conda)"
            echo "[INFO] Conda base: $(conda info --base)"

            if ! conda env list | awk '{{print $1}}' | grep -qx "$ENV_NAME"; then
            echo "[ERROR] Conda env not found: $ENV_NAME"
            conda env list || true
            exit 1
            fi

            cd "$RUN_DIR"
            echo "[INFO] Run directory: $(pwd)"
            echo "[INFO] Executing via conda run: ./run.sh . '$INPUT_FILE' '$OUTPUT_FILE'"

            conda run -n "$ENV_NAME" python -V
            conda run -n "$ENV_NAME" which python || true

            if ! conda run -n "$ENV_NAME" bash ./run.sh . "$INPUT_FILE" "$OUTPUT_FILE"; then
            echo "[ERROR] run.sh failed"
            echo "==================== STDERR (tail) ===================="
            tail -n 200 "$log_err" || true
            echo "======================================================="
            exit 1
            fi

            if [[ -f "$OUTPUT_FILE" ]]; then
            echo "[OK] SUCCESS: output file found: $OUTPUT_FILE"
            exit 0
            fi

            echo "[ERROR] Expected output file not found: $OUTPUT_FILE"
            echo "==================== STDERR (tail) ===================="
            tail -n 200 "$log_err" || true
            echo "==================== STDOUT (tail) ===================="
            tail -n 200 "$log_out" || true
            exit 1
            """
            with open(temp_script_path, "w") as script_file:
                script_file.write(bash_script)

            self.logger.info(f"\nRunning bash script: {temp_script_path}\n")
            with open(temp_script_path, "r") as script_file:
                self.logger.debug(f"Bash script content:\n{script_file.read()}\n")
            try:
                out = run_command(["bash", temp_script_path])
                self.logger.info(f"Bash script output: {out}")
                self.logger.info("Reading output path")
                with open(bash_output_path, "r") as f:
                    self.logger.info(f.read())
                self.logger.info("Done reading output path")
                if os.path.exists(output_log_path):
                    with open(output_log_path, "r") as f:
                        data = f.read()
                        echo(data, fg="cyan", bold=True)
                self.logger.info(f"Bash script subprocess output: {out}")
                logs = read_logs(error_log_path)
                formatted_error = "".join(logs)
                if formatted_error:
                    echo_exceptions(
                        f"Error detected originated from the bash execution: {formatted_error}",
                        ClickInterface(),
                        bg=None,
                        fg="red",
                    )
                bsh_data, bsh_error = read_csv(bash_output_path)
                if bsh_error:
                    raise Exception(f"Could not read bash output: {bsh_error}")
                self.logger.info("Running model for bash data consistency checking")
                if not os.path.exists(input_file_path):
                    raise Exception(
                        "Input file path {0} does not exist".format(
                            os.path.abspath(input_file_path)
                        )
                    )
                cmd = f"ersilia serve {self.model_id} --disable-cache && ersilia -v run -i {os.path.abspath(input_file_path)} -o {output_path}"
                self.logger.info(f"Running command: {cmd}")
                out = run_command(cmd)
                ers_data, ers_error = read_csv(output_path, flag=True)
                if ers_error:
                    raise Exception(f"Could not read ersilia output: {ers_error}")

            except Exception as e:
                return [
                    (
                        (
                            Checks.RUN_BASH.value,
                            f"Detailed error: {e}",
                            str(STATUS_CONFIGS.FAILED),
                        )
                    )
                ]
            self.checkup_service.original_data_list = (
                self.checkup_service._get_original_data_list("csv", input_file_path)
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
            status = _compare_outputs(
                bsh_data,
                ers_data,
                lambda a, b: self._compare_string_similarity(a, b, 95),
            )
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
        metadata = self.ios_service._read_metadata()
        if not str1 or not str2:
            if "Source" in metadata:
                if metadata["Source"] == "Online":
                    return True
        similarity = fuzz.ratio(str1, str2)
        return similarity >= threshold

    def run(self):
        """
        Run the model tests and checks.
        """
        results = []

        def _process_stage(name, method, echo_prefix=True):
            if echo_prefix:
                echo(f"Performing {name} checks.", fg="yellow", bold=True)
            out = method()

            if isinstance(out, tuple) and len(out) == 2:
                good, _bad = out
                results.extend(good)
                sys.exit(1)

            if isinstance(out, list):
                results.extend(out)
            else:
                results.append(out)

        try:
            if not any((self.inspect, self.surface, self.shallow, self.deep)):
                echo(
                    "No flag is specified please at least specify [--inspect].",
                    fg="red",
                    bold=True,
                )
                sys.exit(1)

            self._configure_environment()
            self.setup_service.get_model()

            if self.inspect:
                _process_stage("basic", self._perform_basic_checks)

            if self.surface:
                for name, method in (
                    ("basic", self._perform_basic_checks),
                    ("surface", self._perform_surface_check),
                ):
                    _process_stage(name, method)

            if self.shallow:
                for name, method in (
                    ("basic", self._perform_basic_checks),
                    ("surface", self._perform_surface_check),
                    ("shallow", self._perform_shallow_checks),
                ):
                    _process_stage(name, method)

            if self.deep:
                for name, method in (
                    ("basic", self._perform_basic_checks),
                    ("surface", self._perform_surface_check),
                    ("shallow", self._perform_shallow_checks),
                ):
                    _process_stage(name, method)

                results.append(self._perform_deep_checks())

            self.ios_service.collect_and_save_json(
                results, self.report_file, self.from_dockerhub, self.deep
            )
            echo("Model tests and checks completed.", fg="green", bold=True)
            echo("Deleting model...", fg="yellow", bold=True)
            self.delete()
            echo("Model successfully deleted", fg="green", bold=True)

        except SystemExit as e:
            echo(
                f"Caught SystemExit({e.code}), returning partial results. "
                "Saving report, deleting model and exiting.",
                fg="yellow",
                bold=True,
            )
            self.ios_service.collect_and_save_json(
                results, self.report_file, self.from_dockerhub, self.deep
            )
            self.delete()
            echo("Model successfully deleted", fg="green", bold=True)
            sys.exit(1)

        except Exception as error:
            tb = traceback.format_exc()
            echo(f"An error occurred: {error}\nTraceback:\n{tb}", fg="red", bold=True)
            echo("Deleting model...", fg="yellow", bold=True)
            self.ios_service.collect_and_save_json(
                results, self.report_file, self.from_dockerhub, self.deep
            )
            self.delete()
            echo("Model successfully deleted", fg="green", bold=True)

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
        dim_check = self.checkup_service.check_dim()

        column_name_check = self.checkup_service.check_model_columns_name()
        results.append(
            self._generate_table_from_check(TableType.MODEL_FILE_CHECKS, [dim_check, column_name_check])
        )
       
        results.append(self._log_directory_sizes())

        docker_check = self._docker_yml_column_name_check()
        if isinstance(docker_check, tuple):
            docker_check = docker_check[0]
            results.append(docker_check)
            echo_exceptions(
                "Dependencies are not pinned properly. System is exiting!",
                ClickInterface(),
            )
            return results, 1
        else:
            results.append(docker_check)
        return results

    def _perform_surface_check(self):
        results = []

        out = self.fetch()
        if out.returncode != 0:
            status = [
                (
                    Checks.FETCH_FAILS.value,
                    "Model not fetched successfully",
                    str(STATUS_CONFIGS.FAILED),
                )
            ]
            results.append(
                self._generate_table_from_check(TableType.FETCH_STATUS_SURFACE, status)
            )
            echo_exceptions(
                "Model was not fetched successfully during start of surface checks. System is exiting before proceeding!",
                ClickInterface(),
            )
            return results, 1

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

        self.fetch_output = self._reusable_fetch_output()
        simple_output = self.checkup_service.check_simple_model_output(
            self.run_model, reuse_output=self.fetch_output
        )
        simple_output_async = self.checkup_service.check_simple_model_async_output(
            self.serve_model
        )
        results.append(
            self._generate_table_from_check(TableType.MODEL_RUN_CHECK, simple_output)
        )
        results.append(
            self._generate_table_from_check(
                TableType.ASNC_MODEL_RUN_CHECK, simple_output_async
            )
        )
        if simple_output[0][-1] == str(STATUS_CONFIGS.FAILED):
            echo_exceptions(
                "Model simple run check has problem. System is exiting before proceeding!",
                ClickInterface(),
            )
            return results, 1
        if simple_output_async[0][-1] == str(STATUS_CONFIGS.FAILED):
            echo_exceptions(
                "Model async simple run check has problem. System is exiting before proceeding!",
                ClickInterface(),
            )
            return results, 1

        return results

    def _perform_shallow_checks(self):
        results = []
        metadata = self.ios_service._read_metadata()
        is_online = metadata.get("Source") == "Online"
        model_output = self.checkup_service.check_model_output_content(
            self.run_example, self.run_model
        )
        results.append(
            self._generate_table_from_check(TableType.MODEL_OUTPUT, model_output)
        )
        validations = []
        robustness_output = self.checkup_service.check_robustness(
            self.run_example, self.run_model
        )
        results.append(
            self._generate_table_from_check(
                TableType.ROBUSTNESS_CHECK, robustness_output
            )
        )
        if "Fixed" in self.ios_service.get_output_consistency() and not is_online:
            res = self._run_single_and_example_input_checks()

            validations.append(
                self._generate_table_from_check(TableType.SHALLOW_CHECK_SUMMARY, res)
            )
            bash_results = self.run_bash()
            if self.permissive and self.from_dockerhub:
                bash_results = [
                    (n, d, str(STATUS_CONFIGS.WARNING)) if s == str(STATUS_CONFIGS.FAILED) else (n, d, s)
                    for n, d, s in bash_results
                ]
            validations.append(
                self._generate_table_from_check(
                    TableType.CONSISTENCY_BASH, bash_results
                )
            )
            results.extend(validations)
            if bash_results[0][-1] == str(STATUS_CONFIGS.FAILED):
                echo_exceptions(
                    "Model output is not consistent. System is exiting before proceeding!",
                    ClickInterface(),
                )
                return results, 1

        else:
            row1 = [
                (
                    Checks.MODEL_CONSISTENCY.value,
                    "Skipped for Online source or Varible output consistency",
                    str(STATUS_CONFIGS.SKIPPED),
                )
            ]
            results.append(
                self._generate_table_from_check(TableType.SHALLOW_CHECK_SUMMARY, row1)
            )

            row2 = [
                (
                    Checks.RUN_BASH.value,
                    "Skipped for Online source or Varible output consistency",
                    str(STATUS_CONFIGS.SKIPPED),
                )
            ]
            results.append(
                self._generate_table_from_check(TableType.CONSISTENCY_BASH, row2)
            )

        return results

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
        results = self._generate_table_from_check(
            TableType.DEPENDECY_COLUMN_CHECK, data
        )
        if docker_check_data[0][-1] == str(STATUS_CONFIGS.FAILED):
            return results, 1
        return results

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
