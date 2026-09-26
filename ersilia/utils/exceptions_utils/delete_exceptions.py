from .exceptions import ErsiliaError

# ruff: noqa: D101, D102


class ModelDeleteError(ErsiliaError):
    def __init__(self, model):
        self.message = "Could not delete model {0}.".format(model)
        self.hints = (
            "Check that the model is available locally with 'ersilia catalog --local'."
        )
        ErsiliaError.__init__(self, self.message, self.hints)
