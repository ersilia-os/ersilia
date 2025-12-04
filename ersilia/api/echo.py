try:
    import emoji
except:
    emoji = None
import json
import os
import shutil
import sys

import click
from rich.console import Console
from rich.text import Text

from ..default import SILENCE_FILE
from ..utils.logging import logger
from ..utils.session import get_session_dir

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


def echo(text, harmonize=True, **styles):
    if getattr(logger, "verbosity", 0) == 1:
        return

    if emoji is not None:
        text = emoji.emojize(text, language="alias")

    if not harmonize:
        return click.echo(click.style(text, **styles))

    color = styles.get("fg") or styles.get("color") or styles.get("style")

    if color in ("yellow", "bright_yellow"):
        icon = "⚠"
    elif color in ("red", "bright_red"):
        icon = "✖"
    else:
        icon = "✓"

    text = f"  {icon}  {text}"

    if sys.stdout.isatty():
        width = shutil.get_terminal_size().columns
        text = click.wrap_text(text, width=width)

    if not styles:
        styles["fg"] = "green"

    return click.echo(click.style(text, **styles))


def spinner(text, func, *args, **kwargs):
    if getattr(logger, "verbosity", 0) == 1:
        return func(*args, **kwargs)
    with console.status(Text(text, style="cyan")):
        try:
            result = func(*args, **kwargs)
            console.print(Text(f"  ✓  {text}", style="green"))
            return result
        except Exception:
            console.print(Text(f"  ✖  {text}", style="red"))
            raise
