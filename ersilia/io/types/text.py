import os
import random
import csv
from ...utils.identifiers.protein import TextIdentifier 
from . import EXAMPLES_FOLDER

EXAMPLES = "text.tsv"

class IO(object):
    def __init__(self):
        self.identifier = TextIdentifier()
        self.example_file = os.path.join(EXAMPLES_FOLDER, EXAMPLES)

    def examples(self, n_samples): #TODO
        pass 

    def parse(self, text): #TODO
        data = {"key": text, "input": text, "text": text}
        return data

    def string_delimiter(self):
        return "."

    def column_delimiter(self):
        return ","