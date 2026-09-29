import asyncio
import functools
import logging
import sys

from ... import EOS
from ...default import (
    CURRENT_LOGGING_FILE,
    DEFAULT_ERSILIA_ERROR_EXIT_CODE,
    ERSILIA_CATALOG_URL,
)


def echo(text: str, **styles):
    from ..echo import echo as _echo

    return _echo(text, **styles)


def user_message_and_hints(error):
    # The text users see: the error's own message and hints, never the
    # "Ersilia exception class: ..." block that str(ErsiliaError) holds.
    from ..exceptions_utils.exceptions import ErsiliaError

    if isinstance(error, ErsiliaError):
        message = getattr(error, "message", None)
        hints = getattr(error, "hints", None)
        if not message:
            text = str(error)
            if "Detailed error:\n" in text:
                text = text.split("Detailed error:\n", 1)[1]
            message, _, hints = text.partition("\n\nHints:\n")
        return message.strip(), (hints or "").strip()
    message = str(error).strip()
    return (message or f"Unexpected error ({type(error).__name__})."), ""


def _is_verbose():
    return logging.getLogger("ersilia").level == logging.DEBUG


def _report(error, exit):
    if _is_verbose():
        text = ":police_car_light::police_car_light::police_car_light: Something went wrong with Ersilia :police_car_light::police_car_light::police_car_light:\n"
        echo(text, blink=False, bold=True, fg="red")
        echo("Error message:\n", fg="red", bold=True)
        echo(str(error), fg="red")
    else:
        message, hints = user_message_and_hints(error)
        echo(message, fg="red")
        if hints:
            echo(hints)
    if _is_verbose():
        text = "If this error message is not helpful, open an issue at:\n"
        text += " - https://github.com/ersilia-os/ersilia\n"
        text += "Or feel free to reach out to us at:\n"
        text += " - hello[at]ersilia.io\n\n"
        text += f"Browse the full model catalog at: {ERSILIA_CATALOG_URL}\n\n"
        text += (
            "If you haven't, try to run your command in verbose mode (-v in the CLI)\n"
        )
        text += " - You will find the console log file in: {0}/{1}".format(
            EOS, CURRENT_LOGGING_FILE
        )
        echo(text, fg="green")
    if exit:
        sys.exit(DEFAULT_ERSILIA_ERROR_EXIT_CODE)
    else:
        raise error


def throw_ersilia_exception(exit=True):
    # Ref: https://stackoverflow.com/a/5929165
    def _throw_ersilia_exception(func):
        if asyncio.iscoroutinefunction(func):
            # Errors happen when the coroutine runs, not when it is created,
            # so they must be caught while awaiting it.
            @functools.wraps(func)
            async def async_inner_function(*args, **kwargs):
                try:
                    return await func(*args, **kwargs)
                except Exception as error:
                    _report(error, exit)

            return async_inner_function

        @functools.wraps(func)
        def inner_function(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as error:
                _report(error, exit)
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
