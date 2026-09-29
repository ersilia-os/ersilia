from contextlib import contextmanager

from ..utils.echo import is_quiet, set_quiet
from ..utils.exceptions_utils.throw_ersilia_exception import (
    is_library_mode,
    set_library_mode,
)
from ..utils.logging import logger


@contextmanager
def library_call(verbose=False):
    """
    Run Ersilia code the way a library should.

    Inside this context, errors are raised to the caller (never printed and
    never ending the process), prompts take their default answer instead of
    waiting for input, and log records go to Ersilia's log files only.
    Nothing is printed unless ``verbose`` is True, in which case the same
    progress lines as the CLI are shown.

    Parameters
    ----------
    verbose : bool, optional
        Show the CLI's progress lines.
    """
    previous_quiet = is_quiet()
    previous_library_mode = is_library_mode()
    set_quiet(not verbose)
    set_library_mode(True)
    logger.set_verbosity(0)
    try:
        yield
    finally:
        set_quiet(previous_quiet)
        set_library_mode(previous_library_mode)
