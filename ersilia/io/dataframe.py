import numpy as np
import csv


class Dataframe(object):
    def __init__(self, keys=None, inputs=None, texts=None, values=None, features=None):
        self.keys = keys
        self.inputs = inputs
        self.texts = texts
        self.values = values
        self.features = features
        self._homogenize()

    def iterrows(self):
        for i in range(len(self.keys)):
            result = {
                "key": self.keys[i],
                "input": self.inputs[i],
                "values": self.values[i],
            }
            yield result

    def _homogenize(self):
        # TODO: Deal with other types of inputs
        self.values = np.array(self.values, dtype=np.float32)

    def from_csv(self, filename):
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
