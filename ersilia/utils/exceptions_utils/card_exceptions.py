import os
from .exceptions import ErsiliaError


class CardErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running card command"
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)
