from .exceptions import ErsiliaError


class ServeErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running serve command"
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)
