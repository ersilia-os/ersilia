import json
import os
import re
import shutil
import subprocess
import sys
from collections import namedtuple

from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

try:
    from inputimeout import TimeoutOccurred, inputimeout
except:
    inputimeout = None
    TimeoutOccurred = None

from ..default import _CONDA_BOOTSTRAP, VERBOSE_FILE
from ..utils.logging import make_temp_dir
from ..utils.session import get_session_dir

console = Console()


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


def _looks_like_conda(cmd):
    if isinstance(cmd, str):
        s = cmd.strip()
        return s.startswith("conda ") or " conda " in f" {s} "
    if isinstance(cmd, (list, tuple)) and cmd:
        return cmd[0] == "conda"
    return False


def run_command(cmd, quiet=True):
    """Run a shell command and return a named tuple with stdout, stderr, and return code.
    Parameters
    ----------
    cmd : str or list The command to run.
    quiet : bool, optional Whether to run the command in quiet mode. Defaults to is_quiet().
    """
    if quiet is None:
        quiet = is_quiet()

    use_bash = _looks_like_conda(cmd)

    if isinstance(cmd, str):
        if use_bash:
            script = _CONDA_BOOTSTRAP + "\n" + cmd
            result = subprocess.run(
                ["bash", "-lc", script],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=os.environ,
            )
        else:
            result = subprocess.run(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=os.environ,
            )
    else:
        if use_bash:
            import shlex

            bash_cmd = " ".join(shlex.quote(str(x)) for x in cmd)
            script = _CONDA_BOOTSTRAP + "\n" + bash_cmd
            result = subprocess.run(
                ["bash", "-lc", script],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=os.environ,
            )
        else:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=os.environ,
            )

    CommandResult = namedtuple("CommandResult", ["returncode", "stdout", "stderr"])
    stdout_str = (result.stdout or "").strip()
    stderr_str = (result.stderr or "").strip()
    output = CommandResult(result.returncode, stdout_str, stderr_str)

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


def print_serve_summary(
    model_id,
    slug,
    url,
    pid,
    srv,
    session_dir,
    apis,
    store_stat,
    enable_cache,
    tracking_enabled,
    tracking_use_case,
):
    """
    Print a rich summary table for a served model.
    """
    table = Table(
        show_header=False,
        box=box.SIMPLE,
        expand=False,
        pad_edge=False,
    )

    table.add_row("Model", f"[bold]{model_id}[/bold] ([dim]{slug}[/dim])")
    table.add_row("URL", f"[cyan]{url}[/cyan]")

    if str(pid) != "-1":
        table.add_row("PID", f"[yellow]{pid}[/yellow]")

    table.add_row("Service", f"[yellow]{srv}[/yellow]")
    table.add_row("Session", f"[yellow]{session_dir}[/yellow]")

    all_apis = apis or []
    if "run" in all_apis:
        all_apis = ["run"] + [a for a in all_apis if a != "run"]
    apis_display = ", ".join(all_apis) if all_apis else "run"
    table.add_row("APIs", f"[cyan]{apis_display}[/cyan]")

    table.add_row("Info", "[cyan]info[/cyan]")

    store_style = "red" if store_stat == "Disabled" else "green"
    table.add_row("Isaura Store", f"[{store_style}]{store_stat}[/{store_style}]")

    cache_text = "Enabled" if enable_cache else "Disabled"
    cache_style = "green" if enable_cache else "red"
    table.add_row("Local cache", f"[{cache_style}]{cache_text}[/{cache_style}]")

    if tracking_enabled:
        tracking_text = f"Enabled ({tracking_use_case})"
        tracking_style = "green"
    else:
        tracking_text = "Disabled"
        tracking_style = "red"
    table.add_row("Tracking", f"[{tracking_style}]{tracking_text}[/{tracking_style}]")

    panel = Panel(
        table,
        title="[bold green]Model served[/bold green]",
        expand=False,
        border_style="green",
    )
    console.print()
    console.print(panel)
    console.print()
