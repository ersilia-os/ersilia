try:
    import emoji
except:
    emoji = None
import click
import os
import json
from ..default import EOS, SILENCE_FILE


class Silencer(object):
    def __init__(self):
        self.silence_file = os.path.join(EOS, SILENCE_FILE)
        if not os.path.exists(self.silence_file):
            self.speak()

    def is_silence(self):
        with open(self.silence_file, "r") as f:
            d = json.load(f)
        return d["silence"]

    def speak(self):
        with open(self.silence_file, "w") as f:
            json.dump({"silence": False}, f, indent=4)

    def silence(self):
        with open(self.silence_file, "w") as f:
            json.dump({"silence": True}, f, indent=4)


def echo(text, **styles):
    silencer = Silencer()
    if silencer.is_silence():
        return
    if emoji is not None:
        text = emoji.emojize(text)
    return click.echo(click.style(text, **styles))
