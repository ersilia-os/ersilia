import os
import csv
import random
import importlib
from ...setup.requirements.compound import (
    ChemblWebResourceClientRequirement,
    RdkitRequirement,
)
from ... import logger
from ..readers.file import TabularFileReader
from . import EXAMPLES_FOLDER


EXAMPLES = "compound.tsv"


class IO(object):
    def __init__(self):
        self.logger = logger
        self.example_file = os.path.join(EXAMPLES_FOLDER, EXAMPLES)
        self.setup()
        self.identifier = importlib.import_module(
            "ersilia.utils.identifiers.compound"
        ).CompoundIdentifier()

    def setup(self):
        self.logger.debug("Checking RDKIT and other requirements")
        RdkitRequirement()
        ChemblWebResourceClientRequirement()

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

    def _parse_text(self, text):
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

    def parse(self, data):
        if type(data) is dict:
            return data
        else:
            return self._parse_text(data)

    def is_input(self, text):
        return self.identifier._is_smiles(text)

    def is_key(self, text):
        return self.identifier._is_inchikey(text)
