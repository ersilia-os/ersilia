from .exceptions import ErsiliaError


class DockerImageNotAvailableError(ErsiliaError):
    def __init__(self, model):
        self.message = (
            "Error occured while trying pull docker image of model {0}".format(model)
        )
        self.hints = "Check that the model image ersiliaos/{0} is actually available in Ersilia's DockerHub.\nIf you are working with ARM64 (e.g. M1/M2 Apple chips, it is possible that pulling went wrong because no image with the ARM64 architecture is available".format(
            model
        )
        super().__init__(self.message, self.hints)
