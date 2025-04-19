from .exceptions import ErsiliaError


class NoAWSCredentialsError(ErsiliaError):
    """
    AWS credentials error
    """

    def __init__(self):
        self.message = "You have no AWS credentials in your system. Therefore, tracking is not allowed"
        self.hints = "The variables AWS_ACCESS_KEY_ID and not AWS_SECRET_ACCESS_KEY should be accessible in your system as environmental variables"
        ErsiliaError.__init__(self, self.message, self.hints)
