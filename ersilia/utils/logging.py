import sys
import os
import json
from loguru import logger
from ..default import EOS, LOGGING_FILE, CURRENT_LOGGING_FILE, VERBOSE_FILE


ROTATION = "10 MB"


class Logger(object):
    def __init__(self):
        self.logger = logger
        self.logger.remove()
        self._console = None
        self._file = None
        self.fmt = "{time:HH:mm:ss} | {level: <8} | {message}"
        self._verbose_file = os.path.join(EOS, VERBOSE_FILE)
        self._log_to_console()
        self._log_to_file()
        self._log_to_current_file()
        self._log_terminal_commands_to_console()

    def _log_to_file(self):
        self._file = self.logger.add(
            os.path.join(EOS, LOGGING_FILE), format=self.fmt, rotation=ROTATION
        )

    def _log_to_current_file(self):
        current_log_file = os.path.join(EOS, CURRENT_LOGGING_FILE)

        if os.path.exists(current_log_file):
            os.remove(current_log_file)
        self._current_file = self.logger.add(
            os.path.join(EOS, CURRENT_LOGGING_FILE), format=self.fmt, rotation=ROTATION
        )

    def _log_to_console(self):
        if self._console is None:
            self._console = self.logger.add(sys.stderr, format=self.fmt)

    def _unlog_from_console(self):
        if self._console is not None:
            # self.logger.remove(self._console)
            self.logger.remove()  # TODO: check why console is not found
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
        if verbose:
            self._log_to_console()
            self._log_terminal_commands_to_console()
        else:
            self._unlog_from_console()
            self._unlog_terminal_commands_from_console()
        self.verbosity = verbose

    def debug(self, text):
        self.logger.debug(text)

    def info(self, text):
        self.logger.info(text)

    def warning(self, text):
        self.logger.warning(text)

    def error(self, text):
        self.logger.error(text)

    def critical(self, text):
        self.logger.critical(text)

    def success(self, text):
        self.logger.success(text)


logger = Logger()
