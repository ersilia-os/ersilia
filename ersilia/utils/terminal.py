from subprocess import DEVNULL, STDOUT, check_call, Popen
import os

def run_command(cmd, quiet):
    if type(cmd) == str:
        if quiet:
            with open(os.devnull, "w") as fp:
                Popen(cmd, stdout=fp, shell=True).wait()
        else:
            Popen(cmd, shell=True).wait()
    else:
        if quiet:
            check_call(cmd, stdout=DEVNULL, stderr=STDOUT)
        else:
            check_call(cmd)
