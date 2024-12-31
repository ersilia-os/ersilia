from .exceptions import ErsiliaError


# ruff: noqa: D101
class CardErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running card command"
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)
