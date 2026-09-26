import json
import os
import shutil
import sys

import click
from rich.console import Console
from rich.live import Live
from rich.spinner import Spinner
from rich.text import Text

from ..default import SILENCE_FILE
from ..utils.logging import logger
from ..utils.session import get_session_dir

# Every line is "  <icon>  <text>": the icon at column 2 and the text at
# column 5. The running spinner uses the same columns.
console = Console()


class Silencer(object):
    """
    A class to manage the silencing of CLI output.

    Attributes
    ----------
    silence_file : str
        Path to the silence file.

    Methods
    -------
    is_silence()
        Checks if the CLI output is silenced.
    speak()
        Enables CLI output.
    silence()
        Disables CLI output.
    """

    def __init__(self):
        self.silence_file = os.path.join(get_session_dir(), SILENCE_FILE)
        if not os.path.exists(self.silence_file):
            self.speak()

    def is_silence(self):
        """
        Checks if the CLI output is silenced.

        Returns
        -------
        bool
            True if the CLI output is silenced, False otherwise.
        """
        with open(self.silence_file, "r") as f:
            d = json.load(f)
        return d["silence"]

    def speak(self):
        """
        Enables CLI output.
        """
        with open(self.silence_file, "w") as f:
            json.dump({"silence": False}, f, indent=4)

    def silence(self):
        """
        Disables CLI output.
        """
        with open(self.silence_file, "w") as f:
            json.dump({"silence": True}, f, indent=4)


# One look per kind of message, shared by every command:
#   success  "  ✓  text"  green
#   info     "  ▪  text"  default colour (also progress and hints)
#   warning  "  ⚠  text"  yellow
#   error    "  ✖  text"  red
# Wrapped and multi-line text hangs under the text, at column 5.
_KINDS = {
    "green": ("✓", "green"),
    "bright_green": ("✓", "green"),
    "yellow": ("⚠", "yellow"),
    "bright_yellow": ("⚠", "yellow"),
    "red": ("✖", "red"),
    "bright_red": ("✖", "red"),
}
_INDENT = " " * 5


def _format(text, icon):
    width = shutil.get_terminal_size().columns if sys.stdout.isatty() else None
    lines = []
    for i, line in enumerate(str(text).split("\n")):
        prefix = f"  {icon}  " if i == 0 else _INDENT
        if width and len(prefix) + len(line) > width:
            lines.append(
                click.wrap_text(
                    line,
                    width=width,
                    initial_indent=prefix,
                    subsequent_indent=_INDENT,
                    preserve_paragraphs=False,
                )
            )
        else:
            lines.append(prefix + line)
    return "\n".join(lines)


def echo(text, harmonize=True, **styles):
    """
    Print a message in the CLI's standard style.

    The kind of message is given by ``fg``: green is a success, yellow a
    warning, red an error, and anything else (or nothing) is info.

    Parameters
    ----------
    text : str
        The message. Newlines start indented continuation lines.
    harmonize : bool, optional
        If False, print the text as is, without icon or indentation.
    **styles
        ``fg`` selects the kind of message; ``err=True`` prints to stderr.
        Other click styles are ignored when harmonizing, so every message
        of a kind looks the same.
    """
    if getattr(logger, "verbosity", 0) == 1:
        return
    err = styles.pop("err", False)
    if not harmonize:
        return click.echo(click.style(text, **styles), err=err)
    color = styles.get("fg") or styles.get("color")
    icon, fg = _KINDS.get(color, ("▪", None))
    return click.echo(click.style(_format(text, icon), fg=fg), err=err)


class _SpinnerLine:
    # A running step drawn exactly like a message line, with the spinner in
    # the icon's place: "  ⠋  text".
    def __init__(self, text):
        self._spinner = Spinner("dots", style="green")
        self._text = str(text)

    def __rich_console__(self, console, options):
        frame = self._spinner.render(console.get_time())
        line = Text.assemble("  ", frame, "  ", self._text)
        line.no_wrap = True
        line.overflow = "ellipsis"
        yield line


def _running(text):
    return Live(
        _SpinnerLine(text),
        console=console,
        transient=True,
        refresh_per_second=12,
    )


def confirm(question, default=False):
    """
    Ask a yes/no question, formatted like the other lines.

    Parameters
    ----------
    question : str
        The question, e.g. "Continue?".
    default : bool, optional
        The answer when the user just presses Enter.

    Returns
    -------
    bool
        The answer.
    """
    if not sys.stdin.isatty():
        # Nobody can answer (e.g. a script): use the default and say so,
        # instead of waiting for input that never comes.
        answer = "Yes" if default else "No"
        echo(f"{question} {answer} (no terminal to ask; using the default).")
        return default
    return click.confirm(f"  ▪  {question}", default=default, prompt_suffix=" ")


def spinner(text, func, *args, done=None, **kwargs):
    """
    Run ``func`` while showing ``text``, then print a success line.

    Parameters
    ----------
    text : str
        What is happening, shown while ``func`` runs.
    func : callable
        The work to do.
    done : str, optional
        The success line printed afterwards. Defaults to ``text``.
    """
    if getattr(logger, "verbosity", 0) == 1:
        return func(*args, **kwargs)
    with _running(text):
        try:
            result = func(*args, **kwargs)
        except Exception:
            console.print(Text(_format(text, "✖"), style="red"))
            raise
    console.print(Text(_format(done or text, "✓"), style="green"))
    return result


async def async_spinner(text, coro, done=None):
    """
    Await ``coro`` while showing ``text``, then print a success line.

    Parameters
    ----------
    text : str
        What is happening, shown while ``coro`` runs.
    coro : coroutine
        The work to await.
    done : str, optional
        The success line printed afterwards. Defaults to ``text``.
    """
    if getattr(logger, "verbosity", 0) == 1:
        return await coro
    with _running(text):
        try:
            result = await coro
        except Exception:
            console.print(Text(_format(text, "✖"), style="red"))
            raise
    console.print(Text(_format(done or text, "✓"), style="green"))
    return result
