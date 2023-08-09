import random
import string
import json
import subprocess
import os


def run_command_check_output(cmd):
    return subprocess.check_output(cmd, shell=True, stderr=open(os.devnull))


class ModelIdentifier(object):
    def __init__(self):
        self.letters = string.ascii_lowercase
        self.numbers = "0123456789"

    def encode(self):
        result_str = "eos"
        result_str += random.choice(self.numbers[1:])
        result_str += "".join(
            random.choice(self.letters + self.numbers) for _ in range(3)
        )
        return result_str

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

    def choice(self):
        while True:
            model_id = self.encode()
            if not self.exists(model_id):
                return model_id


if __name__ == "__main__":
    mi = ModelIdentifier()
    print(mi.choice())
