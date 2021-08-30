import random
import string
from ..paths import Paths


class ModelIdentifier(object):
    def __init__(self):
        super().__init__()
        self.letters = string.ascii_lowercase
        self.numbers = "0123456789"

    def encode(self):
        result_str = "eos"
        result_str += random.choice(self.numbers[1:])
        result_str += "".join(
            random.choice(self.letters + self.numbers) for _ in range(3)
        )
        return result_str

    def is_valid(self, s):
        return Paths._eos_regex().match(s)

    def is_test(self, s):
        if s[3] == "0":
            return True
        else:
            return False

    def generate(self, n):
        ids = set()
        while True:
            if len(ids) < n:
                ids.update([self.encode()])
            else:
                break
        ids = list(ids)
        random.shuffle(ids)
        return ids
