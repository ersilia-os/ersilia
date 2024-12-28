import os

from ...utils.identifiers.protein import ProteinIdentifier
from . import EXAMPLES_FOLDER

EXAMPLES = "protein.tsv"

# ruff: noqa: D101, D102


class IO(object):
    def __init__(self):
        self.identifier = ProteinIdentifier()
        self.example_file = os.path.join(EXAMPLES_FOLDER, EXAMPLES)

    def example(self, n_samples):  # TODO
        pass

    def parse(self, text):  # TODO
        pass

    def string_delimiter(self):
        return "."
