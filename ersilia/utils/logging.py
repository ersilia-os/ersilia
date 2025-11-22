import json
import logging
import os
import tempfile
from logging.handlers import RotatingFileHandler
from pathlib import Path

from rich.logging import RichHandler

from ..default import CURRENT_LOGGING_FILE, LOGGING_FILE, VERBOSE_FILE
from ..utils.session import get_session_dir

ROTATION = "10 MB"

# ruff: noqa: D101, D102, F811


def _parse_rotation(rotation_str):
    """
    Parse a rotation string like '10 MB' into bytes.
    """
    try:
        number, unit = rotation_str.split()
        number = float(number)
        unit = unit.upper()
        factor = {
            "B": 1,
            "KB": 1024,
            "MB": 1024 * 1024,
            "GB": 1024 * 1024 * 1024,
        }.get(unit, 1)
        return int(number * factor)
    except Exception:
        # Fallback: no rotation
        return 0


SUCCESS_LEVEL_NUM = 25
logging.addLevelName(SUCCESS_LEVEL_NUM, "SUCCESS")


def make_temp_dir(prefix):
    """
    Create a temporary directory with a symbolic link in the session logs directory.

    Parameters
    ----------
    prefix : str
        The prefix for the temporary directory name.

    Returns
    -------
    str
        The path to the temporary directory.
    """
    tmp_dir = tempfile.mkdtemp(prefix=prefix)
    tmp_dirname = os.path.basename(tmp_dir)
    logs_tmp_dir = os.path.join(get_session_dir(), "logs", "tmp")
    if not os.path.exists(logs_tmp_dir):
        os.makedirs(logs_tmp_dir)
    dst = Path(os.path.join(logs_tmp_dir, tmp_dirname))
    src = Path(tmp_dir)
    dst.symlink_to(src, target_is_directory=True)
    return tmp_dir


class Logger(object):
    """
    A singleton class to manage logging in Ersilia.

    Methods
    -------
    set_verbosity(verbose)
        Set the verbosity of the logger.
    debug(text)
        Log a debug message.
    info(text)
        Log an info message.
    warning(text)
        Log a warning message.
    error(text)
        Log an error message.
    critical(text)
        Log a critical message.
    success(text)
        Log a success message.
    """

    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(Logger, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        self.logger = logging.getLogger("ersilia")
        self.logger.setLevel(logging.INFO)
        # Clear any existing handlers to mimic loguru's .remove()
        self.logger.handlers.clear()

        self._console = None
        self._file = None
        self._current_file = None
        self.fmt = "%(asctime)s | %(levelname)-8s | %(message)s"

        self._verbose_file = os.path.join(get_session_dir(), VERBOSE_FILE)
        self._log_to_console()
        self._log_to_file()
        self._log_to_current_file()
        self._log_terminal_commands_to_console()

    def _log_to_file(self):
        logging_file = os.path.join(get_session_dir(), LOGGING_FILE)
        if self._file is None and os.path.exists(logging_file):
            max_bytes = _parse_rotation(ROTATION)
            if max_bytes > 0:
                handler = RotatingFileHandler(
                    logging_file, maxBytes=max_bytes, backupCount=5, mode="w"
                )
            else:
                handler = logging.FileHandler(logging_file, mode="w")
            handler.setFormatter(logging.Formatter(self.fmt))
            self.logger.addHandler(handler)
            self._file = handler

    def _log_to_current_file(self):
        session_dir = get_session_dir()
        current_log_file = os.path.join(session_dir, CURRENT_LOGGING_FILE)
        if self._current_file is None:
            Path(current_log_file).parent.mkdir(parents=True, exist_ok=True)
            max_bytes = _parse_rotation(ROTATION)
            try:
                if max_bytes > 0:
                    handler = RotatingFileHandler(
                        current_log_file, maxBytes=max_bytes, backupCount=5
                    )
                else:
                    handler = logging.FileHandler(current_log_file)
            except OSError:
                return
            handler.setFormatter(logging.Formatter(self.fmt))
            self.logger.addHandler(handler)
            self._current_file = handler

    def _log_to_console(self):
        if self._console is None:
            rich_handler = RichHandler(
                rich_tracebacks=True,
                markup=True,
                log_time_format="%H:%M:%S",
                show_path=False,
            )
            formatter = logging.Formatter("%(message)s")
            rich_handler.setFormatter(formatter)
            self.logger.addHandler(rich_handler)
            self._console = rich_handler

    def _unlog_from_console(self):
        if self._console is not None:
            try:
                self.logger.removeHandler(self._console)
            except Exception:
                pass
            self._console = None

    def _log_terminal_commands_to_console(self):
        d = {"verbose": True}
        with open(self._verbose_file, "w") as f:
            json.dump(d, f, indent=4)

    def _unlog_terminal_commands_from_console(self):
        d = {"verbose": False}
        with open(self._verbose_file, "w") as f:
            json.dump(d, f, indent=4)

    def set_verbosity(self, verbose):
        """
        Set the verbosity of the logger.

        Parameters
        ----------
        verbose : bool
            Whether to enable verbose logging.
        """
        if verbose:
            self.logger.setLevel(logging.DEBUG)  # <-- add this
            self._log_to_console()
            self._log_terminal_commands_to_console()
        else:
            self.logger.setLevel(logging.INFO)  # <-- and this
            self._unlog_from_console()
            self._unlog_terminal_commands_from_console()
        self.verbosity = verbose

    def debug(self, text):
        """
        Log a debug message.

        Parameters
        ----------
        text : str
            The message to log.
        """
        self.logger.debug(text)

    def info(self, text):
        """
        Log an info message.

        Parameters
        ----------
        text : str
            The message to log.
        """
        self.logger.info(text)

    def warning(self, text):
        """
        Log a warning message.

        Parameters
        ----------
        text : str
            The message to log.
        """
        self.logger.warning(text)

    def error(self, text):
        """
        Log an error message.

        Parameters
        ----------
        text : str
            The message to log.
        """
        self.logger.error(text)

    def critical(self, text):
        """
        Log a critical message.

        Parameters
        ----------
        text : str
            The message to log.
        """
        self.logger.critical(text)

    def success(self, text):
        """
        Log a success message.

        Parameters
        ----------
        text : str
            The message to log.
        """
        self.logger.log(SUCCESS_LEVEL_NUM, text)


logger = Logger()
