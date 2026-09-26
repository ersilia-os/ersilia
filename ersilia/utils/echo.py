import json
import os
import shutil
import sys

import click
from rich.console import Console
from rich.live import Live
from rich.padding import Padding
from rich.panel import Panel
from rich.progress import BarColumn, Progress, TextColumn, TimeElapsedColumn
from rich.spinner import Spinner
from rich.table import Table
from rich.text import Text
from rich.theme import Theme

from ..default import SILENCE_FILE
from ..utils.logging import logger
from ..utils.session import get_session_dir

# The look of everything the CLI prints, in one place.
#
# Layout: every line is "  <icon>  <text>": the icon at column 2 and the text
# at column 5. Continuation lines, progress bars and the running spinner use
# the same columns, and panels are indented by 2 so their border lines up with
# the icons.
#
# Colour means something, so it is used sparingly: green for success, yellow
# for warnings, red for errors, cyan for links. Everything else is the
# terminal's default colour, with secondary details dimmed. Nothing is bold.
THEME = Theme(
    {
        "bar.back": "bright_black",
        "bar.complete": "default",
        "bar.finished": "green",
        "bar.pulse": "default",
        "progress.elapsed": "dim",
        "progress.remaining": "dim",
        "progress.download": "dim",
        "progress.percentage": "dim",
        "progress.spinner": "default",
        "status.spinner": "default",
    }
)
console = Console(theme=THEME, highlight=False)


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
        self._spinner = Spinner("dots")
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


def progress_bar(*columns):
    """
    A progress bar that lines up under the text of the line above it.

    Parameters
    ----------
    *columns : rich.progress.ProgressColumn
        Columns shown after the bar (e.g. a count), before the elapsed time.

    Returns
    -------
    rich.progress.Progress
        Use it as a context manager, like any rich Progress.
    """
    return Progress(
        TextColumn("    "),
        BarColumn(bar_width=40),
        *columns,
        TimeElapsedColumn(),
        console=console,
    )


def fields_table():
    """
    A two-column table of fields and values, for panels.

    Returns
    -------
    rich.table.Table
        Labels are dimmed; add rows with ``table.add_row(label, value)``.
    """
    table = Table(show_header=False, box=None, padding=(0, 2), pad_edge=False)
    table.add_column(style="dim", no_wrap=True)
    table.add_column(overflow="fold")
    return table


def print_panel(renderable, title=None):
    """
    Print a panel in the CLI's style: indented, with a quiet border.

    Parameters
    ----------
    renderable : rich renderable
        The panel's content, usually a ``fields_table()``.
    title : str, optional
        Plain text shown in the top border.
    """
    panel = Panel(
        renderable,
        title=Text(title) if title else None,
        title_align="left",
        border_style="bright_black",
        expand=False,
        padding=(1, 2),
    )
    console.print()
    console.print(Padding(panel, (0, 0, 0, 2), expand=False))
    console.print()


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
