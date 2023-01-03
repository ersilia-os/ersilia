from .exceptions import ErsiliaError


class ModelDeleteError(ErsiliaError):
    def __init__(self, model):
        self.message = "Error occured while deleting model {0}".format(model)
        self.hints = "Check that the model is actually installed in your local device:\n$ ersilia serve {0}".format(
            model
        )
        super().__init__(self.message, self.hints)
