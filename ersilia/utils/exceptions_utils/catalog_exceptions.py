from .exceptions import ErsiliaError


class CatalogErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running catalog command"
        self.hints = ""
        super().__init__(self.message, self.hints)
