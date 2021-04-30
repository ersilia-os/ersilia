import random
import os
from ...utils.identifiers.protein import ProteinIdentifier
from . import EXAMPLES_FOLDER

EXAMPLES = "protein.tsv"


class ProteinIO(object):
    def __init__(self):
        self.identifier = ProteinIdentifier()
        self.example_file = os.path.join(EXAMPLES_FOLDER, EXAMPLES)

    def random(self):
        with open(self.example_file, "r") as f:
            line = random.choice(f.readlines()).rstrip("\n").split("\t")
        prot = {"uniprot_ac": line[0], "seq": line[1], "name": line[2]}
        return prot

    def single(self, inp):  # TODO
        pass

    def multiple(self, filename):  # TODO
        pass
