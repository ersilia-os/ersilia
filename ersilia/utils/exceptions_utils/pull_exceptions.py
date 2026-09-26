from .exceptions import ErsiliaError

# ruff: noqa: D101, D102


class DockerImageNotAvailableError(ErsiliaError):
    def __init__(self, model):
        self.message = "Could not pull the Docker image of model {0}.".format(model)
        self.hints = "Check that the image ersiliaos/{0} exists on DockerHub. On ARM64 machines (e.g. Apple M1/M2), the model may have no ARM64 image.".format(
            model
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class DockerImageArchitectureNotAvailableError(ErsiliaError):
    def __init__(self, model):
        self.message = "Could not pull model {0} from DockerHub.".format(model)
        self.hints = "The model may not support your machine's architecture (e.g. Apple M1/M2). You can run Ersilia in GitHub Codespaces instead, or contact us at hello@ersilia.io for help."
        ErsiliaError.__init__(self, self.message, self.hints)


class DockerConventionalPullError(ErsiliaError):
    def __init__(self, model):
        self.message = "Could not pull the Docker image of model {0}.".format(model)
        ErsiliaError.__init__(self, self.message)
