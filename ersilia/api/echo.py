try:
    import emoji
except:
    emoji = None
import json
import os

import click

from ..default import SILENCE_FILE
from ..utils.session import get_session_dir


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


def echo(text, **styles):
    """
    Prints styled text to the CLI, optionally with emojis, unless silenced.

    Parameters
    ----------
    text : str
        The text to be printed.
    **styles : dict
        Additional styling options for the text.

    Returns
    -------
    None
    """
    silencer = Silencer()
    if silencer.is_silence():
        return
    if emoji is not None:
        text = emoji.emojize(text)
    return click.echo(click.style(text, **styles))
