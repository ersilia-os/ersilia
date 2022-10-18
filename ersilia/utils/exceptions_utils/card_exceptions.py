from .exceptions import ErsiliaError

class CardErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running card command"
        self.hints = ""
        super().__init__(self.message, self.hints)