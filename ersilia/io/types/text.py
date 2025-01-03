import csv
import importlib
import os
import random

from ... import logger
from ...utils.identifiers.arbitrary import ArbitraryIdentifier
from ..shape import InputShapeList, InputShapePairOfLists, InputShapeSingle
from . import EXAMPLES_FOLDER
from .examples import text as test_examples

EXAMPLES = "text.tsv"


class IO(object):
    """
    Class to handle input/output operations for text data.

    Parameters
    ----------
    input_shape : object
        Input shape specification.
    """

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

    def setup(self):
        """
        Setup necessary requirements for text inputs.

        Returns
        -------
        None
        """
        self.logger.debug("No setup is necessary for text file")

    def example(self, n_samples):
        """
        Generate example data.

        Parameters
        ----------
        n_samples : int
            Number of samples to generate.

        Returns
        -------
        generator
            Generator yielding example data.
        """
        return self._example(n_samples)

    def test(self):
        """
        Get test examples.

        Returns
        -------
        list
            List of test examples.
        """
        return self._test

    def parse(self, datum):
        """
        Parse a datum.

        Parameters
        ----------
        datum : any
            The datum to parse.

        Returns
        -------
        dict
            Parsed datum.
        """
        if type(datum) is dict:
            return datum
        else:
            return self._parser(datum)

    def is_input(self, text):
        """
        Check if a text is a valid input.

        Parameters
        ----------
        text : str
            The text to check.

        Returns
        -------
        bool
            True if the text is a valid input, False otherwise.
        """
        if text == "input":
            return False
        if self.identifier._is_checksum(text):
            return False
        return True  # TODO

    def is_key(self, text):
        """
        Check if a text is a valid key.

        Parameters
        ----------
        text : str
            The text to check.

        Returns
        -------
        bool
            True if the text is a valid key, False otherwise.
        """
        if text == "key":
            return False
        return self.identifier._is_checksum(text)

    def string_delimiter(self):
        """
        Get the string delimiter.

        Returns
        -------
        str
            The string delimiter.
        """
        # TODO
        return "$$$"

    def column_delimiter(self):
        """
        Get the column delimiter.

        Returns
        -------
        str
            The column delimiter.
        """
        # TODO
        return ","

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
