from .exceptions import ErsiliaError

# ruff: noqa: D101, D102


class ServeErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running serve command"
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)


class BadGatewayError(ErsiliaError):
    def __init__(self, url):
        self.message = "The url {0} returned a 502 (bad gateway) error".format(url)
        self.hints = "If you are trying to access a docker container from this url ({0}), check that the architecture matches your system's architecture. If this is not the case, it is possible that the model simply cannot run on your system.".format(
            url
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class DockerNotActiveError(ErsiliaError):
    def __init__(self):
        self.message = "Docker is not active. Cannot serve model"
        self.hints = "Please activate docker and try again"
        ErsiliaError.__init__(self, self.message, self.hints)
