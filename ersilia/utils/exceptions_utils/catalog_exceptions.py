from .exceptions import ErsiliaError


class CatalogErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running catalog command"
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)
