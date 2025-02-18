class BentoMLException(Exception):
    """
    Exception raised for errors in the BentoML tool.

    Parameters
    ----------
    message : str
        A custom error message describing the issue.
    """

    def __init__(self, message: str):
        super().__init__(message)


class BentoMLConfigException(Exception):
    """
    Exception raised for configuration errors in the BentoML tool.
    """

    pass
