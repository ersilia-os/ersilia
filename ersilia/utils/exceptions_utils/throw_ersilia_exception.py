import asyncio
import functools
import logging
import sys

from ...default import (
    DEFAULT_ERSILIA_ERROR_EXIT_CODE,
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


# In library mode (the Python API) errors are re-raised for the caller to
# handle: nothing is printed and the process never exits.
_library_mode = False


def set_library_mode(enabled):
    """
    Make decorated functions re-raise errors instead of printing and exiting.

    Parameters
    ----------
    enabled : bool
        True for library (Python API) behaviour.
    """
    global _library_mode
    _library_mode = bool(enabled)


def is_library_mode():
    """
    Tell whether errors are re-raised for the caller.

    Returns
    -------
    bool
        True in library mode.
    """
    return _library_mode


def _is_verbose():
    return logging.getLogger("ersilia").level == logging.DEBUG


class ErsiliaErrorReported(Exception):
    """
    An error that was already shown to the user.

    Raised instead of ``sys.exit`` inside coroutines: a ``SystemExit`` raised
    in an asyncio task is also kept on the task, and asyncio later prints it
    with a traceback ("Task exception was never retrieved"). Whoever runs the
    coroutine turns this into the exit code, without printing anything more.

    Parameters
    ----------
    code : int
        The exit code.
    """

    def __init__(self, code=DEFAULT_ERSILIA_ERROR_EXIT_CODE):
        self.code = code
        Exception.__init__(
            self, "Ersilia error already reported (exit {0})".format(code)
        )


def show_error(error):
    """
    Show an error the standard way: the message, then its hints.

    Always shown, also in verbose mode (where it follows the log lines).

    Parameters
    ----------
    error : Exception
        The error.
    """
    message, hints = user_message_and_hints(error)
    echo(message, fg="red", force=True)
    if hints:
        echo(hints, force=True)
    if _is_verbose():
        echo(
            "If this does not help, please open an issue at https://github.com/ersilia-os/ersilia/issues",
            force=True,
        )


def _report(error, exit, in_coroutine=False):
    show_error(error)
    if exit:
        if in_coroutine:
            raise ErsiliaErrorReported() from None
        sys.exit(DEFAULT_ERSILIA_ERROR_EXIT_CODE)
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
                except ErsiliaErrorReported:
                    raise
                except Exception as error:
                    if _library_mode:
                        raise
                    _report(error, exit, in_coroutine=True)

            return async_inner_function

        @functools.wraps(func)
        def inner_function(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except ErsiliaErrorReported as reported:
                # Shown by a coroutine this function ran.
                if exit:
                    sys.exit(reported.code)
                raise
            except Exception as error:
                if _library_mode:
                    raise
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
