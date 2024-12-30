import hashlib


class TextIdentifier(object):
    """
    A class to handle text identification by generating MD5 checksums.

    This class provides methods to generate a unique identifier (checksum) for a given text string using the MD5 hashing algorithm.
    It also includes a method to perform a basic validation check on the generated checksum.
    """

    def __init__(self):
        pass

    @staticmethod
    def _is_checksum(text):
        # TODO this method is not secure, for now, it is just a quick check
        if not text.startswith("key"):
            return False
        if len(text) != 32 + 3:
            return False
        text = text[3:]
        if " " in text:
            return False
        if "," in text:
            return False
        if not all(c in "0123456789abcdefABCDEF" for c in text):
            return False
        return True

    def encode(self, text: str) -> str:
        """
        Generate an MD5 checksum for the given text.

        This method takes a text string as input and returns a unique identifier (checksum) for the text using the MD5 hashing algorithm.
        The checksum is prefixed with "key" to distinguish it from other strings.

        Parameters
        ----------
        text : str
            The text string to generate a checksum for.

        Returns
        -------
        str
            The MD5 checksum of the text, prefixed with "key".
        """
        return "key" + hashlib.md5(text.encode("utf-8")).hexdigest()


Identifier = TextIdentifier
