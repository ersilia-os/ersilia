import subprocess
import os


def run_command(cmd, quiet):
    if type(cmd) == str:
        if quiet:
            with open(os.devnull, "w") as fp:
                subprocess.Popen(cmd, stdout=fp, shell=True, env=os.environ).wait()
        else:
            subprocess.Popen(cmd, shell=True, env=os.environ).wait()
    else:
        if quiet:
            subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, env=os.environ)
        else:
            subprocess.check_call(cmd, env=os.environ)


def run_command_check_output(cmd):
    result = subprocess.run(cmd, stdout=subprocess.PIPE, env=os.environ)
    return result.stdout
