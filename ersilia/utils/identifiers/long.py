try:
    import uuid
except ModuleNotFoundError as err:
    uuid = None
import random

ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
PATTERN = [8, 4, 4, 4, 12]
SEP = "-"


class LongIdentifier(object):
    def __init__(self):
        pass

    @staticmethod
    def encode():
        """Get UUID code (long identifier)"""
        if uuid is None:
            alphabet = ALPHABET.lower()
            for n in PATTERN:
                s += ["".join([random.choice(alphabet) for _ in range(n)])]
            return "-".join(s)
        else:
            return str(uuid.uuid4())
