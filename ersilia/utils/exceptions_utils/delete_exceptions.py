from .exceptions import ErsiliaError

class DeleteErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running delete command"
        self.hints = ""
        super().__init__(self.message, self.hints)