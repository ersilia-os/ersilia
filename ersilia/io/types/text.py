import os
import csv
import random
import importlib

from ...utils.identifiers.arbitrary import ArbitraryIdentifier
from ... import logger
from ..shape import InputShapeSingle, InputShapeList, InputShapePairOfLists
from .examples import text as test_examples
from . import EXAMPLES_FOLDER


EXAMPLES = "text.tsv"


class IO(object):
    def __init__(self, input_shape):
        self.logger = logger
        self.input_shape = input_shape
        self.example_file = os.path.join(EXAMPLES_FOLDER, EXAMPLES)
        self.setup()
        self.identifier = importlib.import_module(
            "ersilia.utils.identifiers.text"
        ).TextIdentifier()
        self.arbitrary_identifier = ArbitraryIdentifier()
        if type(self.input_shape) is InputShapeSingle:
            self.logger.debug(
                "InputShapeSingle shape: {0}".format(self.input_shape.name)
            )
            self._example = self._example_single
            self._parser = self._parse_single
            self._test = test_examples.input_shape_single_text

        if type(self.input_shape) is InputShapeList:
            # TODO
            self.logger.warning(
                "Input shape list for text is not available in Ersilia yet!"
            )
            raise Exception
            self.logger.debug("InputShapeList shape: {0}".format(self.input_shape.name))
            self._example = self._example_list
            self._parser = self._parse_list
            self._test = test_examples.input_shape_list_text

        if type(self.input_shape) is InputShapePairOfLists:
            # TODO
            self.logger.warning(
                "Input shape pair of list for text is not available in Ersilia yet!"
            )
            raise Exception

    def _sample_example_singlets(self, n_samples):
        delimiter = None
        if self.example_file.endswith(".tsv"):
            delimiter = "\t"
        if self.example_file.endswith(".csv"):
            delimiter = ","
        with open(self.example_file, "r") as f:
            reader = csv.reader(f, delimiter=delimiter)
            R = []
            for r in reader:
                R += [r]
        idxs = [i for i in range(len(R))]
        idxs = random.choices(idxs, k=n_samples)
        D = []
        for idx in idxs:
            r = R[idx]
            D += [{"key": r[0], "input": r[2], "text": r[2]}]
        return D

    def _example_single(self, n_samples):
        D = self._sample_example_singlets(n_samples)
        for d in D:
            yield d

    def _example_list(self, n_samples):
        D = []
        for _ in range(n_samples):
            D_ = self._sample_example_singlets(10)
            input = [x["input"] for x in D_]
            text = self.string_delimiter().join(input)
            key = self.arbitrary_identifier.encode(text)
            D += [{"key": key, "input": input, "text": text}]
        for d in D:
            yield d

    def setup(self):
        self.logger.debug("No setup is necessary for text file")

    def example(self, n_samples):
        return self._example(n_samples)

    def test(self):
        return self._test

    def _parse_text(self, datum):
        text = str(datum)
        inp = text
        key = self.identifier.encode(inp)
        result = {"key": key, "input": inp, "text": text}
        return result

    def _parse_single(self, datum):
        return self._parse_text(datum)

    def _parse_list(self, datum):
        inp = datum
        text = self.string_delimiter().join(inp)
        key = self.arbitrary_identifier.encode(text)
        result = {"key": key, "input": inp, "text": text}
        return result

    def parse(self, datum):
        if type(datum) is dict:
            return datum
        else:
            return self._parser(datum)

    def is_input(self, text):
        if text == "input":
            return False
        if self.identifier._is_checksum(text):
            return False
        return True  # TODO

    def is_key(self, text):
        if text == "key":
            return False
        return self.identifier._is_checksum(text)

    def string_delimiter(self):
        # TODO
        return "$$$"

    def column_delimiter(self):
        # TODO
        return ","
