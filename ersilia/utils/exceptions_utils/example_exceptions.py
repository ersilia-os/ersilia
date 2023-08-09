from .exceptions import ErsiliaError


class ExampleErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running example command"
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)
