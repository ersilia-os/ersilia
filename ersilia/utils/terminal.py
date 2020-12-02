import subprocess
import os


def run_command(cmd, quiet):
    if type(cmd) == str:
        if quiet:
            with open(os.devnull, "w") as fp:
                subprocess.Popen(cmd, stdout=fp, shell=True).wait()
        else:
            subprocess.Popen(cmd, shell=True).wait()
    else:
        if quiet:
            subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        else:
            subprocess.check_call(cmd)


def run_command_check_output(cmd):
    result = subprocess.run(cmd, stdout=subprocess.PIPE)
    return result.stdout
