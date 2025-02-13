import csv
import json
import os
import shutil
import subprocess
import sys

try:
    from inputimeout import TimeoutOccurred, inputimeout
except:
    inputimeout = None
    TimeoutOccurred = None

from collections import namedtuple

from ..default import OUTPUT_DATASTRUCTURE, VERBOSE_FILE
from ..utils.logging import make_temp_dir
from ..utils.session import get_session_dir
from .hdf5 import Hdf5DataLoader


def is_quiet():
    """
    Check if the current session is in quiet mode.

    Returns
    -------
    bool
        True if the session is in quiet mode, False otherwise.
    """
    verbose_file = os.path.join(get_session_dir(), VERBOSE_FILE)
    if not os.path.exists(verbose_file):
        return False
    with open(verbose_file, "r") as f:
        d = json.load(f)
    if d["verbose"]:
        return False
    else:
        return True


def run_command(cmd, quiet=None):
    """
    Run a shell command and return a named tuple with stdout, stderr, and return code.

    Parameters
    ----------
    cmd : str or list
        The command to run.
    quiet : bool, optional
        Whether to run the command in quiet mode. Defaults to is_quiet().
    """

    if quiet is None:
        quiet = is_quiet()

    # Run the command and capture outputs
    result = subprocess.run(
        cmd,
        shell=isinstance(cmd, str),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=os.environ,
    )

    CommandResult = namedtuple("CommandResult", ["returncode", "stdout", "stderr"])
    stdout_str = result.stdout.strip()
    stderr_str = result.stderr.strip()
    output = CommandResult(
        returncode=result.returncode, stdout=stdout_str, stderr=stderr_str
    )

    # Log outputs if not in quiet mode
    if not quiet:
        if stdout_str:
            print(stdout_str)
        if stderr_str:
            print(stderr_str, file=sys.stderr)

    return output


def run_command_check_output(cmd):
    """
    Run a shell command and capture its output.

    Parameters
    ----------
    cmd : str or list
        The command to run.

    Returns
    -------
    str
        The output of the command.
    """
    if type(cmd) is str:
        assert ">" not in cmd
        tmp_folder = make_temp_dir(prefix="ersilia-")
        os.listdir(tmp_folder)
        tmp_file = os.path.join(tmp_folder, "out.txt")
        cmd = cmd + " > " + tmp_file
        with open(os.devnull, "w") as fp:
            subprocess.Popen(cmd, shell=True, stderr=fp, stdout=fp).wait()
        with open(tmp_file, "r") as f:
            result = f.read()
        shutil.rmtree(tmp_folder)
        return result
    else:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, env=os.environ)
        return result.stdout


def raw_input_with_timeout(prompt, default_answer, timeout=5):
    """
    Prompt the user for input with a timeout.

    Parameters
    ----------
    prompt : str
        The prompt message.
    default_answer : str
        The default answer if the user does not respond within the timeout.
    timeout : int, optional
        The timeout in seconds. Default is 5.

    Returns
    -------
    str
        The user's input or the default answer if the timeout is reached.
    """
    if inputimeout is None:
        return input(prompt)
    try:
        answer = inputimeout(prompt=prompt, timeout=timeout)
    except TimeoutOccurred:
        answer = default_answer
    return answer


def yes_no_input(prompt, default_answer, timeout=5):
    """
    Prompt the user for a yes/no input with a timeout.

    Parameters
    ----------
    prompt : str
        The prompt message.
    default_answer : str
        The default answer if the user does not respond within the timeout.
    timeout : int, optional
        The timeout in seconds. Default is 5.

    Returns
    -------
    bool
        True if the user's input is 'yes', False otherwise.
    """
    ans = raw_input_with_timeout(
        prompt=prompt, default_answer=default_answer, timeout=timeout
    )
    if ans is None or ans == "":
        ans = default_answer
    ans = str(ans).lower()
    if ans[0] == "n":
        return False
    else:
        return True


def _flatten_data(json_data):
    flattened = []
    for item in json_data:
        row = {}
        for key, value in item.items():
            if isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    if sub_key == "text":
                        continue
                    row[sub_key] = _handle_supported_structures(sub_value)
            elif key != "text":
                row[key] = _handle_supported_structures(value)
        flattened.append(row)
    return flattened


def _handle_supported_structures(value):
    for dtype, checker in OUTPUT_DATASTRUCTURE.items():
        if checker(value):
            if dtype == "Single":
                return str(value[0])
            elif dtype == "List" or dtype == "Flexible List":
                return ", ".join(map(str, value))
            elif dtype == "Matrix":
                return "; ".join(", ".join(map(str, row)) for row in value)
            elif dtype == "Serializable Object":
                return json.dumps(value)
    return str(value)


def _read_hdf5_with_loader(file_path):
    loader = Hdf5DataLoader()
    loader.load(file_path)

    data = []
    for i, key in enumerate(loader.keys):
        row = {
            "Key": key,
            "Input": loader.inputs[i] if i < len(loader.inputs) else None,
            "Value": loader.values[i] if i < len(loader.values) else None,
            "Feature": loader.features[i] if i < len(loader.features) else None,
        }
        data.append(row)
    return data


def _read_csv(file_path):
    with open(file_path, mode="r") as file:
        reader = csv.DictReader(file)
        return [dict(row) for row in reader]


def print_result_table(data):
    """
    Print a result table with solid borders from JSON, CSV, or HDF5-like data.
    Supports formatted JSON strings.
    """
    COLOR, BORDER_CHAR, VERTICAL_BORDER = "\033[0m", "━", "┃"

    if isinstance(data, str):
        try:
            parsed_data = json.loads(data)
            if isinstance(parsed_data, list) and all(
                isinstance(item, dict) for item in parsed_data
            ):
                data = parsed_data
            else:
                raise ValueError(
                    "The JSON string must represent a list of dictionaries."
                )
        except json.JSONDecodeError:
            if os.path.isfile(data):
                if data.endswith(".json"):
                    with open(data, mode="r") as file:
                        data = json.load(file)
                elif data.endswith(".csv"):
                    data = _read_csv(data)
                elif data.endswith((".h5", ".hdf5")):
                    data = _read_hdf5_with_loader(data)
                else:
                    raise ValueError(f"Unsupported file type: {data}")
            else:
                raise ValueError(
                    f"Provided string is neither valid JSON nor a file path. {data}"
                )

    if isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict):
        if data[0].get("input") or data[0].get("output"):
            data = _flatten_data(data)

        headers = list(data[0].keys())

        column_widths = {
            header: max(len(header), max(len(str(row.get(header, ""))) for row in data))
            + 2
            for header in headers
        }

        def format_row(row_data, is_header=False):
            if is_header:
                return (
                    COLOR
                    + VERTICAL_BORDER
                    + VERTICAL_BORDER.join(
                        f" {header.ljust(column_widths[header])} " for header in headers
                    )
                    + VERTICAL_BORDER
                    + COLOR
                )
            else:
                return (
                    COLOR
                    + VERTICAL_BORDER
                    + VERTICAL_BORDER.join(
                        f" {str(row_data.get(header, '')).ljust(column_widths[header])} "
                        for header in headers
                    )
                    + VERTICAL_BORDER
                    + COLOR
                )

        total_width = sum(column_widths.values()) + (3 * len(headers))
        border_line = BORDER_CHAR * total_width

        print(border_line)
        print(format_row(headers, is_header=True))
        print(border_line)
        for row in data:
            print(format_row(row))
            print(border_line)
    else:
        print(
            f"Invalid input data format. Please provide valid JSON, CSV, or HDF5 data.{data}"
        )
