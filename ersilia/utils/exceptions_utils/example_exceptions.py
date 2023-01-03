from .exceptions import ErsiliaError


class ExampleErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running example command"
        self.hints = ""
        super().__init__(self.message, self.hints)
