try:
    import uuid
except ModuleNotFoundError:
    uuid = None
import random

ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
PATTERN = [8, 4, 4, 4, 12]
SEP = "-"


class LongIdentifier(object):
    """
    A class to generate long identifiers (UUIDs).

    Methods
    -------
    encode()
        Generate a UUID or a random identifier if UUID is not available.
    """

    def __init__(self):
        pass

    @staticmethod
    def encode():
        """
        Generate a UUID or a random identifier if UUID is not available.

        Returns
        -------
        str
            A UUID string or a randomly generated identifier.
        """
        if uuid is None:
            alphabet = ALPHABET.lower()
            s = []
            for n in PATTERN:
                s += ["".join([random.choice(alphabet) for _ in range(n)])]
            return "-".join(s)
        else:
            return str(uuid.uuid4())


Identifier = LongIdentifier
