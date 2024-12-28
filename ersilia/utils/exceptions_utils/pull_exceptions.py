from .exceptions import ErsiliaError

# ruff: noqa: D101, D102


class DockerImageNotAvailableError(ErsiliaError):
    def __init__(self, model):
        self.message = (
            "Error occured while trying pull docker image of model {0}".format(model)
        )
        self.hints = "Check that the model image ersiliaos/{0} is actually available in Ersilia's DockerHub.\nIf you are working with ARM64 (e.g. M1/M2 Apple chips, it is possible that pulling went wrong because no image with the ARM64 architecture is available".format(
            model
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class DockerImageArchitectureNotAvailableError(ErsiliaError):
    def __init__(self, model):
        self.message = "It was not possible to pull model {0} from Ersilia's DockerHub repository.".format(
            model
        )
        self.hints = "If you are using an Apple M1/M2 chip, it is possible that this model is not supported for your architecture, unfortunately.\nOne possible alternative is to use GitHub Codespaces to run Ersilia on the cloud, and fetch the model from there. If you absolutely want this model to run on a Mac, please reach out to us and we will try to help."
        ErsiliaError.__init__(self, self.message, self.hints)


class DockerConventionalPullError(ErsiliaError):
    def __init__(self, model):
        self.message = "Conventional pull did not work for model {0}".format(model)
        ErsiliaError.__init__(self, self.message)
