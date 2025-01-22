import csv

import numpy as np


class Dataframe(object):
    """
    Class to handle data in a dataframe-like structure.

    Parameters
    ----------
    keys : list, optional
        List of keys.
    inputs : list, optional
        List of inputs.
    texts : list, optional
        List of texts.
    values : list, optional
        List of values.
    features : list, optional
        List of features.
    """

    def __init__(self, keys=None, inputs=None, texts=None, values=None, features=None):
        self.keys = keys
        self.inputs = inputs
        self.texts = texts
        self.values = values
        self.features = features
        self._homogenize()

    def iterrows(self):
        """
        Iterate over the rows of the dataframe.

        Yields
        ------
        dict
            Dictionary representing a row in the dataframe.
        """
        for i in range(len(self.keys)):
            result = {
                "key": self._process(self.keys, i),
                "input": self._process(self.inputs, i),
                "text": self._process(self.texts, i),
                "values": self._process(self.values, i),
            }
            yield result

    def from_csv(self, filename):
        """
        Load data from a CSV file.

        Parameters
        ----------
        filename : str
            Path to the CSV file.

        Returns
        -------
        None
        """
        keys = []
        inputs = []
        values = []
        with open(filename, "r") as f:
            reader = csv.reader(f)
            h = next(reader)
            for r in reader:
                keys += [r[0]]
                inputs += [r[1]]
                values += [r[2:]]
            features = h[2:]
        self.keys = keys
        self.inputs = inputs
        self.texts = None
        self.values = values
        self.features = features
        self._homogenize()

    def _process(self, col, idx):
        if col is None:
            return None
        else:
            return col[idx]

    def _float(self, x):
        try:
            return float(x)
        except:
            return np.nan

    def _homogenize(self):
        if self.values is None:
            return
        values = np.zeros((len(self.values), len(self.values[0])))
        for i, v in enumerate(self.values):
            for j, x in enumerate(v):
                values[i, j] = self._float(x)
        self.values = np.array(values, dtype=np.float32)
        # TODO: Deal with other types of inputs
