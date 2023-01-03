from .exceptions import ErsiliaError


class CloseErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running close command"
        self.hints = ""
        super().__init__(self.message, self.hints)
