import h5py
import numpy as np


class Hdf5Data(object):
    """
    A class to handle HDF5 data storage.

    Parameters
    ----------
    values : array-like
        The data values.
    keys : array-like
        The keys associated with the data.
    inputs : array-like
        The inputs associated with the data.
    features : array-like
        The features associated with the data.
    """
    def __init__(self, values, keys, inputs, features):
        self.values = np.array(values, dtype=np.float32)
        self.keys = np.array(keys, dtype=h5py.string_dtype())
        self.inputs = np.array(inputs, dtype=h5py.string_dtype())
        self.features = np.array(features, dtype=h5py.string_dtype())

    def save(self, filename):
        """
        Save the data to an HDF5 file.

        Parameters
        ----------
        filename : str
            The path to the HDF5 file.
        """
        with h5py.File(filename, "w") as f:
            f.create_dataset("Values", data=self.values)
            f.create_dataset("Keys", data=self.keys)
            f.create_dataset("Inputs", data=self.inputs)
            f.create_dataset("Features", data=self.features)


class Hdf5DataLoader(object):
    """
    A class to load data from HDF5 files.

    Methods
    -------
    load(h5_file)
        Load data from an HDF5 file.
    """
    def __init__(self):
        self.values = None
        self.keys = None
        self.inputs = None
        self.features = None

    def load(self, h5_file):
        """
        Load data from an HDF5 file.

        Parameters
        ----------
        h5_file : str
            The path to the HDF5 file.
        """
        with h5py.File(h5_file, "r") as f:
            self.values = f["Values"][:]
            self.keys = [x.decode("utf-8") for x in f["Keys"][:]]
            self.inputs = [x.decode("utf-8") for x in f["Inputs"][:]]
            self.features = [x.decode("utf-8") for x in f["Features"][:]]


class Hdf5DataStacker(object):
    """
    A class to stack multiple HDF5 files into one.

    Parameters
    ----------
    h5_files : list
        A list of paths to the HDF5 files to stack.
    """
    def __init__(self, h5_files):
        self.h5_files = h5_files

    def stack(self, h5_file):
        """
        Stack the HDF5 files into one.

        Parameters
        ----------
        h5_file : str
            The path to the output HDF5 file.
        """
        with h5py.File(h5_file, "a") as f0:
            for i, h5 in enumerate(self.h5_files):
                with h5py.File(h5, "r") as f1:
                    values = f1["Values"][:]
                    keys = f1["Keys"][:]
                    inputs = f1["Inputs"][:]
                    features = f1["Features"][:]
                    if i == 0:
                        f0.create_dataset(
                            "Values",
                            data=values,
                            shape=values.shape,
                            maxshape=(None, values.shape[1]),
                            chunks=True,
                        )
                        f0.create_dataset(
                            "Keys",
                            data=keys,
                            shape=keys.shape,
                            maxshape=(None,),
                            chunks=True,
                        )
                        f0.create_dataset(
                            "Inputs",
                            data=inputs,
                            shape=inputs.shape,
                            maxshape=(None,),
                            chunks=True,
                        )
                        f0.create_dataset(
                            "Features",
                            data=features,
                            shape=features.shape,
                            maxshape=(None,),
                            chunks=True,
                        )
                    else:
                        f0["Values"].resize(
                            (f0["Values"].shape[0] + values.shape[0]), axis=0
                        )
                        f0["Values"][-values.shape[0] :] = values
                        f0["Keys"].resize((f0["Keys"].shape[0] + keys.shape[0]), axis=0)
                        f0["Keys"][-keys.shape[0] :] = keys
                        f0["Inputs"].resize(
                            (f0["Inputs"].shape[0] + inputs.shape[0]), axis=0
                        )
                        f0["Inputs"][-inputs.shape[0] :] = inputs
