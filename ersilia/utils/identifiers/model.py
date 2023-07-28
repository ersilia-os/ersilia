import random
import string
import json
from ..paths import Paths
from ..terminal import run_command_check_output


class ModelIdentifier(object):
    def __init__(self):
        pass
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
        if len(s) != 7:
            return False
        return Paths._eos_regex().match(s)

    def is_test(self, s):
        if s[3] == "0":
            return True
        else:
            return False

    def exists(self, model_id):
        cmd = "curl https://api.github.com/repos/ersilia-os/{0}".format(model_id)
        output = json.loads(run_command_check_output(cmd))
        if "message" in output:
            if output["message"] == "Not Found":
                return False
            else:
                return True
        else:
            return False

    def generate(self, n):
        ids = set()
        while True:
            if len(ids) < n:
                model_id = self.encode()
                if not self.exists(model_id):
                    ids.update([self.encode()])
            else:
                break
        ids = list(ids)
        random.shuffle(ids)
        return ids

    def choice(self):
        while True:
            model_id = self.encode()
            if not self.exists(model_id):
                return model_id
