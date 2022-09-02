
import sys

from ...cli import echo
from email_reporting import send_exception_report_email

def query_yes_no(question):
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    prompt = " [y/n] "

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' " "(or 'y' or 'n').\n")

try:
    raise Exception # this is where the code for handling the CLI input should go
except Exception as E:
 
    text = ":triangular_flag: Something went wrong with Ersilia...\n\n"
    text += "{}\n\n".format(self.__class__.__name__)
    echo(text)
    echo("Error message:\n")
    echo(":prohibited: " + str(E), fg="red")
    text = "If this error message is not helpful, open an issue at:\n"
    text += " - https://github.com/ersilia-os/ersilia\n"
    text += "Or feel free to reach out to us at:\n"
    text += " - hello[at]ersilia.io\n\n"
    text += "If you haven't, try to run your command in verbose mode (-v in the CLI)\n\n"
    echo(text)
 
    if query_yes_no("Would you like to report this error to Ersilia?"):
        send_exception_report_email(E)

    if query_yes_no("Would you like to access the log?"):
        # execute cli logic for [y/n] query and write log to a file