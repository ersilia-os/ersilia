from .exceptions import ErsiliaError


class CloseErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running close command"
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)
