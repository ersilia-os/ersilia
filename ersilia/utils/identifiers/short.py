import random

try:
    from hashids import Hashids
except ModuleNotFoundError:
    Hashids = None
try:
    from datetime import datetime
except ModuleNotFoundError:
    datetime = None

ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
RANDINT = 9999999
LENGTH = 8


class ShortIdentifier(object):
    """
    A class to generate short identifiers.
    """

    def __init__(self):
        if Hashids is None:
            self.hashids = None
        else:
            self.hashids = Hashids(salt="ersilia is open source", alphabet=ALPHABET)

    def encode(self):
        """
        Generate a short identifier based on the current timestamp or a random number.

        Returns
        -------
        str
            A short identifier string.
        """
        if self.hashids is None:
            return "".join([random.choice(ALPHABET) for _ in range(LENGTH)])
        else:
            if datetime is None:
                number = random.randint(RANDINT)
            else:
                number = int(datetime.today().timestamp())
            return str(self.hashids.encode(number))


Identifier = ShortIdentifier
