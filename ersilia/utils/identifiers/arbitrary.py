import hashlib


class ArbitraryIdentifier(object):
    """
    Class for handling arbitrary identifiers.
    """

    def __init__(self):
        pass

    def encode(self, text: str) -> str:
        """
        Encode the given text using MD5.

        Parameters
        ----------
        text : str
            The text to encode.

        Returns
        -------
        str
            The encoded text.
        """
        return hashlib.md5(text.encode("utf-8")).hexdigest()


Identifier = ArbitraryIdentifier
