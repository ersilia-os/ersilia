import random

try:
    from hashids import Hashids
except ModuleNotFoundError as err:
    Hashids = None
try:
    from datetime import datetime
except ModuleNotFoundError as err:
    datetime = None

ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
RANDINT = 9999999
LENGTH = 8


class ShortIdentifier(object):
    def __init__(self):
        if Hashids is None:
            self.hashids = None
        else:
            self.hashids = Hashids(salt="ersilia is open source", alphabet=ALPHABET)

    def encode(self):
        """Short identifier based on timestamp"""
        if self.hashids is None:
            return "".join([random.choice(ALPHABET) for _ in range(LENGTH)])
        else:
            if datetime is None:
                number = random.randint(RANDINT)
            else:
                number = int(datetime.today().timestamp())
            return str(self.hashids.encode(number))
