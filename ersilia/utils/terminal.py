import os
import subprocess
import json
from ..default import EOS, VERBOSE_FILE


def is_quiet():
    verbose_file = os.path.join(EOS, VERBOSE_FILE)
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
    result = subprocess.run(cmd, stdout=subprocess.PIPE, env=os.environ)
    return result.stdout
