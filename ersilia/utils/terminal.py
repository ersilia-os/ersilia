import json
import os
import re
import psutil
import shutil
import subprocess
import sys
import hashlib
import time
import uuid
from pathlib import Path

try:
    from inputimeout import TimeoutOccurred, inputimeout
except:
    inputimeout = None
    TimeoutOccurred = None

from collections import namedtuple

from ..utils.logging import make_temp_dir
from ..utils.session import get_session_dir

from ..default import VERBOSE_FILE, STATE_DIRECTORY


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


def is_quoted_list(s: str) -> bool:
    pattern = r"^(['\"])\[.*\]\1$"
    return bool(re.match(pattern, s))


def get_terminal_session_id(state_dir=STATE_DIRECTORY):
    try:
        tty = os.ttyname(0)
    except Exception:
        tty = "notty"

    try:
        sid = os.getsid(0)
    except Exception:
        sid = os.getsid(os.getpid())

    key_plain = f"{tty}|sid:{sid}"
    key = hashlib.sha1(key_plain.encode()).hexdigest()  # compact stable key

    # 2) Persist/lookup the ID
    if state_dir is None:
        # XDG-ish default, works on macOS/Linux. Falls back to ~/.local/state
        state_root = os.environ.get("XDG_STATE_HOME") or str(Path.home() / ".local" / "state")
        state_dir = os.path.join(state_root, "terminal-session-ids")
    os.makedirs(state_dir, exist_ok=True)

    index_path = os.path.join(state_dir, "index.json")
    try:
        with open(index_path, "r", encoding="utf-8") as f:
            index = json.load(f)
    except FileNotFoundError:
        index = {}

    entry = index.get(key)
    if entry:
        return entry["id"]

    # New session => mint a new ID and record minimal metadata
    term_id = str(uuid.uuid4())
    index[key] = {
        "id": term_id,
        "tty": tty,
        "sid": sid,
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    # Best-effort write (atomic-ish)
    tmp = index_path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(index, f, indent=2)
    os.replace(tmp, index_path)

    return term_id


def is_terminal_session_alive(term_id, state_dir=STATE_DIRECTORY):
    
    if state_dir is None:
        state_root = os.environ.get("XDG_STATE_HOME") or str(Path.home() / ".local" / "state")
        state_dir = os.path.join(state_root, "terminal-session-ids")

    index_path = os.path.join(state_dir, "index.json")
    try:
        with open(index_path, "r", encoding="utf-8") as f:
            index = json.load(f)
    except FileNotFoundError:
        return False

    entry = None
    for key, val in index.items():
        if val["id"] == term_id:
            entry = val
            break

    if not entry:
        return False

    tty, sid = entry["tty"], entry["sid"]

    if not os.path.exists(tty):
        return False

    for proc in psutil.process_iter(attrs=["pid"]):
        try:
            if os.getsid(proc.info["pid"]) == sid:
                return True
        except Exception:
            continue
    return False
