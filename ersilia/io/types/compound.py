import random
import os
from ...utils.identifiers.compound import CompoundIdentifier
from ...setup.requirements.rdkit import RdkitRequirement
from . import EXAMPLES_FOLDER
from ... import logger

EXAMPLES = "compound.tsv"


class IO(object):
    def __init__(self):
        self.logger = logger
        self.identifier = CompoundIdentifier()
        self.example_file = os.path.join(EXAMPLES_FOLDER, EXAMPLES)
        self.setup()

    def setup(self):
        self.logger.debug("Checking RDKIT requirement")
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

    def is_input(self, text):
        return self.identifier._is_smiles(text)

    def is_key(self, text):
        return self.identifier._is_inchikey(text)
