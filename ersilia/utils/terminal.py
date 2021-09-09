import os
import subprocess
import json
from ..default import EOS, SILENCE_FILE


def is_quiet():
    silence_file = os.path.join(EOS, SILENCE_FILE)
    with open(silence_file, "r") as f:
        d = json.load(f)
    return d["silence"]


def run_command(cmd):
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
    result = subprocess.run(cmd, stdout=subprocess.PIPE, env=os.environ)
    return result.stdout
