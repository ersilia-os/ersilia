import sys

import click

from ... import EOS
from ...default import CURRENT_LOGGING_FILE, DEFAULT_ERSILIA_ERROR_EXIT_CODE

try:
    import emoji
except:
    emoji = None


def echo(text: str, **styles):
    if emoji is not None:
        text = emoji.emojize(text)
    return click.echo(click.style(text, **styles))


def throw_ersilia_exception(exit=True):
    # Ref: https://stackoverflow.com/a/5929165
    def _throw_ersilia_exception(func):
        def inner_function(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as error:
                text = ":police_car_light::police_car_light::police_car_light: Something went wrong with Ersilia :police_car_light::police_car_light::police_car_light:\n"
                echo(text, blink=False, bold=True, fg="red")
                echo("Error message:\n", fg="red", bold=True)
                echo(str(error), fg="red")
                text = "If this error message is not helpful, open an issue at:\n"
                text += " - https://github.com/ersilia-os/ersilia\n"
                text += "Or feel free to reach out to us at:\n"
                text += " - hello[at]ersilia.io\n\n"
                text += "If you haven't, try to run your command in verbose mode (-v in the CLI)\n"
                text += " - You will find the console log file in: {0}/{1}".format(
                    EOS, CURRENT_LOGGING_FILE
                )
                echo(text, fg="green")
                if exit:
                    sys.exit(DEFAULT_ERSILIA_ERROR_EXIT_CODE)
                else:
                    raise error
                # FIXME: Enable automatic reporting of issues
                # if query_yes_no("Would you like to report this error to Ersilia?"):
                #     if query_yes_no(
                #         "Would you like to include your last Ersilia command in the issue (for issue reproducibility)?"
                #     ):
                #         sys.stdout.write("Please re-type your last Ersilia command: ")
                #         message = input()
                #         send_exception_issue(E, message)
                #     else:
                #         send_exception_issue(E, "")

                # if query_yes_no("Would you like to access the log?"):
                #     print("No log info")
                #     # TODO: execute cli logic for [y/n] query and write log to a file

        return inner_function

    return _throw_ersilia_exception
