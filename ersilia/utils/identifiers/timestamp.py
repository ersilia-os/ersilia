from datetime import datetime


class TimeStampIdentifier(object):
    """
    Class for handling timestamp identifiers.
    """

    def __init__(self):
        self.stamp = datetime.now()

    def encode(self):
        """
        Encode the current timestamp.

        Returns
        -------
        str
            The encoded timestamp.
        """
        return self.stamp.strftime("%Y%m%d%H%M%S")


Identifier = TimeStampIdentifier
