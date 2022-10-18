from .exceptions import ErsiliaError


class ServeErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running serve command"
        self.hints = ""
        super().__init__(self.message, self.hints)
