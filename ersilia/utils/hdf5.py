import h5py
import numpy as np


class Hdf5Data(object):
    def __init__(self, values, keys, inputs, features):
        self.values = np.array(values, dtype=np.float32)
        self.keys = np.array(keys, dtype=h5py.string_dtype())
        self.inputs = np.array(inputs, dtype=h5py.string_dtype())
        self.features = np.array(features, dtype=h5py.string_dtype())

    def save(self, filename):
        with h5py.File(filename, "w") as f:
            f.create_dataset("Values", data=self.values)
            f.create_dataset("Keys", data=self.keys)
            f.create_dataset("Inputs", data=self.inputs)
            f.create_dataset("Features", data=self.features)


class Hdf5DataLoader(object):
    def __init__(self):
        self.values = None
        self.keys = None
        self.inputs = None
        self.features = None

    def load(self, h5_file):
        with h5py.File(h5_file, "r") as f:
            self.values = f["Values"][:]
            self.keys = [x.decode("utf-8") for x in f["Keys"][:]]
            self.inputs = [x.decode("utf-8") for x in f["Inputs"][:]]
            self.features = [x.decode("utf-8") for x in f["Features"][:]]
