import random
import os
from ...utils.identifiers.compound import CompoundIdentifier
from ...setup.requirements.rdkit import RdkitRequirement
from . import EXAMPLES_FOLDER

EXAMPLES = "compound.tsv"


class IO(object):
    def __init__(self):
        self.identifier = CompoundIdentifier()
        self.example_file = os.path.join(EXAMPLES_FOLDER, EXAMPLES)
        self.setup()

    def setup(self):
        RdkitRequirement()

    def random(self):
        with open(self.example_file, "r") as f:
            line = random.choice(f.readlines()).rstrip("\n").split("\t")
        result = {"key": line[0], "input": line[1], "text": line[2]}
        return result

    def parse(self, text):
        text_type = self.identifier.guess_type(text)
        key = None
        if text_type == "smiles":
            inp = text
        elif text_type == "inchikey":
            inp = self.identifier.unichem_resolver(text)
            key = text
        else:
            inp = self.identifier.chemical_identifier_resolver(text)
        if key is None:
            key = self.identifier.encode(inp)
        result = {"key": key, "input": inp, "text": text}
        return result
