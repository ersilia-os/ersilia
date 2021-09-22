import os
import csv
import random
from ...utils.identifiers.compound import CompoundIdentifier
from ...setup.requirements.rdkit import RdkitRequirement
from ... import logger
from ..readers.file import TabularFileReader
from . import EXAMPLES_FOLDER


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

    def example(self, n_samples):
        tfr = TabularFileReader(self)
        tfr.get_delimiter(self.example_file)
        R = []
        with open(self.example_file, "r") as f:
            reader = csv.reader(f, delimiter=tfr.delimiter)
            for r in reader:
                R += [r]
        idxs = [i for i in range(len(R))]
        idxs = random.choices(idxs, k=n_samples)
        for idx in idxs:
            r = R[idx]
            d = {"key": r[0], "input": r[1], "text": r[2]}
            yield d

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
