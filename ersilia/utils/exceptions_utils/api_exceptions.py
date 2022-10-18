from .exceptions import ErsiliaError

class ApiErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running api command"
        self.hints = ""
        super().__init__(self.message, self.hints)