from .exceptions import ErsiliaError

# ruff: noqa: D101, D102


class ModelDeleteError(ErsiliaError):
    def __init__(self, model):
        self.message = "Error occured while deleting model {0}".format(model)
        self.hints = "Check that the model is actually installed in your local device:\n$ ersilia serve {0}".format(
            model
        )
        ErsiliaError.__init__(self, self.message, self.hints)
