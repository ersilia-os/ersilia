from .exceptions import ErsiliaError

# ruff: noqa: D101, D102


class CatalogErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running catalog command"
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)
