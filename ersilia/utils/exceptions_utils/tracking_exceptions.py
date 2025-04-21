from .exceptions import ErsiliaError


class NoAWSCredentialsError(ErsiliaError):
    """
    AWS credentials error
    """

    def __init__(self):
        self.message = "You have no valid AWS credentials in your system. Therefore, tracking is not allowed"
        self.hints = "The variables AWS_ACCESS_KEY_ID and not AWS_SECRET_ACCESS_KEY should be accessible in your system as environmental variables or in the ~/.aws/credentials file"
        ErsiliaError.__init__(self, self.message, self.hints)


class TrackingNotSupportedError(ErsiliaError):
    """
    Tracking not supported error
    """

    def __init__(self):
        self.message = "Tracking is not supported for this model as currently fetch, or with the specified input and output"
        self.hints = "Tracking is only supported for models that are served from a pulled Docker image, and for which input and outputs are specified as CSV files"
        ErsiliaError.__init__(self, self.message, self.hints)
