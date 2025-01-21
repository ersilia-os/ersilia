import csv
import io
import json
import os
import shutil
import subprocess

from .logging import logger

try:
    from inputimeout import TimeoutOccurred, inputimeout
except:
    inputimeout = None
    TimeoutOccurred = None

from ..default import VERBOSE_FILE
from ..utils.logging import make_temp_dir
from ..utils.session import get_session_dir


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
    Run a shell command.

    Parameters
    ----------
    cmd : str or list
        The command to run.
    quiet : bool, optional
        Whether to run the command in quiet mode. Default is None.
    """
    if quiet is None:
        quiet = is_quiet()
    if type(cmd) == str:
        if quiet:
            with open(os.devnull, "w") as fp:
                subprocess.Popen(
                    cmd, stdout=fp, stderr=fp, shell=True, env=os.environ
                ).wait()
        else:
            subprocess.Popen(cmd, shell=True, env=os.environ).wait()
    else:
        if quiet:
            subprocess.check_call(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                env=os.environ,
            )
        else:
            subprocess.check_call(cmd, env=os.environ)


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


def truncate_output(output, max_length=10):
    """
    Truncate the given output for better readability in the terminal.

    Handles strings, lists, dictionaries, integers, floats, CSV, HDF5 files, and matrices.

    Parameters
    ----------
    output : str, list, dict, int, float, or file object
        The output to truncate.
    max_length : int, optional
        The maximum number of characters or items to display. Default is 10.

    Returns
    -------
    str
        The truncated output as a string.
    """
    if isinstance(output, str):
        return (
            output if len(output) <= 100 else output[:100] + "..."
        )  # Strings: 100 chars
    elif isinstance(output, (int, float)):
        return str(output)  # Numbers: Show as-is
    elif isinstance(output, list):
        if all(isinstance(row, list) for row in output):  # Matrix
            truncated_rows = [
                row[:max_length] for row in output[:max_length]
            ]  # Truncate rows/columns
            return json.dumps(truncated_rows, indent=2) + (
                "\n..." if len(output) > max_length else ""
            )
        else:  # Regular list
            return json.dumps(output[:max_length], indent=2) + (
                "..." if len(output) > max_length else ""
            )
    elif isinstance(output, dict):  # Serializable Object
        truncated = list(output.items())[:max_length]
        return json.dumps(truncated, indent=2) + (
            "..." if len(output) > max_length else ""
        )
    elif hasattr(output, "read"):  # Handle file-like objects
        if output.name.endswith(".csv"):
            output.seek(0)
            reader = csv.reader(output)
            rows = [row for _, row in zip(range(max_length), reader)]
            return "\n".join([",".join(row) for row in rows]) + (
                "\n..." if len(list(reader)) > max_length else ""
            )
        elif output.name.endswith((".h5", ".hdf5")):
            import h5py

            with h5py.File(output.name, "r") as f:
                keys = list(f.keys())[:max_length]
                return f"HDF5 file: {output.name}\nDatasets: {keys}" + (
                    "\n..." if len(keys) > max_length else ""
                )
    else:
        return str(output)


def print_result_table(data):
    """
    Print a result table from CSV or JSON-like data.

    Parameters
    ----------
    data : str or list
        The path to a CSV file or JSON-like data.
    """
    HEADER_COLOR = "\033[95m"
    ROW_COLOR = "\033[94m"
    RESET_COLOR = "\033[0m"
    if isinstance(data, str) and os.path.isfile(data):
        with open(data, mode="r") as file:
            reader = csv.DictReader(file)
            data = [dict(row) for row in reader]
    elif isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict):
        headers = list(data[0].keys())
        column_widths = {
            header: max(len(header), max(len(str(row[header])) for row in data)) + 5
            for header in headers
        }

        def format_row(row_data, is_header=False):
            if is_header:
                return (
                    HEADER_COLOR
                    + " | ".join(
                        f"{header.ljust(column_widths[header])}" for header in headers
                    )
                    + RESET_COLOR
                )
            else:
                return (
                    ROW_COLOR
                    + " | ".join(
                        f"{str(row_data[header]).ljust(column_widths[header])}"
                        for header in headers
                    )
                    + RESET_COLOR
                )

        separator = "-" * (sum(column_widths.values()) + (3 * len(headers) - 1))
        print(separator)
        print(format_row(headers, is_header=True))
        print(separator)
        for row in data:
            print(format_row(row))
        print(separator)
    else:
        logger.debug(
            "Invalid input data format. Please provide either a CSV file path or JSON-like data."
        )


def read_csv_from_string(csv_string):
    """
    Read CSV data from a string.

    Parameters
    ----------
    csv_string : str
        The CSV data as a string.

    Returns
    -------
    list
        A list of dictionaries representing the CSV data.
    """
    f = io.StringIO(csv_string)
    reader = csv.DictReader(f)
    return [dict(row) for row in reader]
