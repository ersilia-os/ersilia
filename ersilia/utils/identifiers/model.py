import json
import random
import string

from ..paths import Paths
from ..terminal import run_command_check_output


class ModelIdentifier(object):
    """
    A class to handle model identification generation for new ersilia model and validation.

    """

    def __init__(self):
        self.letters = string.ascii_lowercase
        self.numbers = "0123456789"

    def encode(self):
        """
        Generate a new model identifier.

        Returns
        -------
        str
            A new model identifier.
        """
        result_str = "eos"
        result_str += random.choice(self.numbers[1:])
        result_str += "".join(
            random.choice(self.letters + self.numbers) for _ in range(3)
        )
        return result_str

    def is_valid(self, s):
        """
        Check if a given string is a valid model identifier.

        Parameters
        ----------
        s : str
            The string to check.

        Returns
        -------
        bool
            True if the string is a valid model identifier, False otherwise.
        """
        if len(s) != 7:
            return False
        return Paths._eos_regex().match(s)

    def is_test(self, s):
        """
        Check if a given model identifier is a test identifier.

        Parameters
        ----------
        s : str
            The model identifier to check.

        Returns
        -------
        bool
            True if the model identifier is a test identifier, False otherwise.
        """
        if s[3] == "0":
            return True
        else:
            return False

    def exists(self, model_id):
        """
        Check if a model identifier exists in the Ersilia repository.

        Parameters
        ----------
        model_id : str
            The model identifier to check.

        Returns
        -------
        bool
            True if the model identifier exists, False otherwise.
        """
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
        """
        Generate a list of unique model identifiers.

        Parameters
        ----------
        n : int
            The number of model identifiers to generate.

        Returns
        -------
        list
            A list of unique model identifiers.
        """
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
        """
        Generate a unique model identifier that does not exist in the Ersilia repository.

        Returns
        -------
        str
            A unique model identifier.
        """
        while True:
            model_id = self.encode()
            if not self.exists(model_id):
                return model_id


Identifier = ModelIdentifier
