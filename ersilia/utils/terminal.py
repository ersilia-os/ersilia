import os
import subprocess
import json
import tempfile
import shutil
import csv
from .logging import logger
import io

try:
    from inputimeout import inputimeout, TimeoutOccurred
except:
    inputimeout = None
    TimeoutOccurred = None

from ..default import VERBOSE_FILE
from ..utils.session import get_session_dir
from ..utils.logging import make_temp_dir


def is_quiet():
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
    if inputimeout is None:
        return input(prompt)
    try:
        answer = inputimeout(prompt=prompt, timeout=timeout)
    except TimeoutOccurred:
        answer = default_answer
    return answer


def yes_no_input(prompt, default_answer, timeout=5):
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
    

def print_result_table(data):
    HEADER_COLOR = "\033[95m"
    ROW_COLOR = "\033[94m"
    RESET_COLOR = "\033[0m"
    if isinstance(data, str) and os.path.isfile(data):
        with open(data, mode="r") as file:
            reader = csv.DictReader(file)
            data = [dict(row) for row in reader]
    if isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict):
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
    f = io.StringIO(csv_string)
    reader = csv.DictReader(f)
    return [dict(row) for row in reader]
