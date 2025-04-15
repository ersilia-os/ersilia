import json

import h5py
import numpy as np


class Hdf5Data:
    """
    A class to handle HDF5 data storage.

    Parameters
    ----------
    values : array-like or nested array-like
        The data values.  Can be a 1D list or a list-of-lists.
    keys : array-like
        The keys associated with the data.
    inputs : array-like
        The inputs associated with the data.
    features : array-like
        The features associated with the data.
    dtype : {float, int, str, None}
        If float/int/str: force all `values` to that type (float→np.float32,
        int→np.int32, str→utf-8 HDF5 string).  If None: infer per your original logic.
    """

    def __init__(self, values, keys, inputs, features, dtype):
        self.dtype = dtype
        if dtype is None:
            self._force_dtype = False
            self._np_dtype = None
        else:
            self._force_dtype = True
            if dtype is float:
                self._np_dtype = np.float32
            elif dtype is int:
                self._np_dtype = np.int32
            elif dtype is str:
                self._np_dtype = h5py.string_dtype(encoding="utf-8")
            else:
                raise ValueError("dtype must be float, int, str, or None")

        self.values = self._convert_values(values, inputs)

        str_dt = h5py.string_dtype(encoding="utf-8")
        self.keys = np.array(keys, dtype=str_dt)
        self.inputs = np.array(inputs, dtype=str_dt)
        self.features = np.array(features, dtype=str_dt)

    def _default_fill(self):
        if np.issubdtype(self._np_dtype, np.floating):
            return np.nan
        elif np.issubdtype(self._np_dtype, np.integer):
            return 0
        else:
            return ""

    def _convert_values(self, values, inputs):
        if not self._force_dtype:
            return self._infer_values(values, inputs)

        if not values:
            fill = self._default_fill()
            return np.full((len(inputs),), fill, dtype=self._np_dtype)

        if all(isinstance(v, (list, tuple)) for v in values):
            rows = [self._convert_1d(row, inputs) for row in values]
            try:
                arr = np.stack(rows)
            except ValueError:
                arr = np.array(rows, dtype=self._np_dtype)
        else:
            arr = self._convert_1d(values, inputs)

        return arr.astype(self._np_dtype, copy=False)

    def _convert_1d(self, values, inputs):
        if not values:
            fill = self._default_fill()
            return np.full((len(inputs),), fill, dtype=self._np_dtype)

        fill = self._default_fill()
        cleaned = [fill if v is None else v for v in values]

        if self._np_dtype == h5py.string_dtype(encoding="utf-8"):
            cleaned = [json.dumps(v) if not isinstance(v, str) else v for v in cleaned]

        return np.array(cleaned, dtype=self._np_dtype)

    def _infer_values(self, values, inputs):
        if not values:
            return np.array(["" for _ in range(len(inputs))], dtype=np.float32)

        if all(isinstance(v, (list, tuple)) for v in values):
            converted_inner = [
                self._convert_values_1d(inner, inputs) for inner in values
            ]
            first_dtype = converted_inner[0].dtype
            if all(arr.dtype == first_dtype for arr in converted_inner):
                try:
                    return np.stack(converted_inner)
                except Exception:
                    return np.array(converted_inner, dtype=first_dtype)
            else:
                converted_inner = [
                    arr.astype(h5py.string_dtype(encoding="utf-8"))
                    for arr in converted_inner
                ]
                return np.array(
                    converted_inner, dtype=h5py.string_dtype(encoding="utf-8")
                )
        else:
            return self._convert_values_1d(values, inputs)

    def _convert_values_1d(self, values, inputs):
        if not values:
            return np.array(
                ["" for _ in range(len(inputs))],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )

        if all(v is None or not v for v in values):
            values = [None for _ in values]
            return np.array(values, dtype=np.float32)

        if all(isinstance(v, str) and v == "" for v in values):
            return np.array(values, dtype=h5py.string_dtype(encoding="utf-8"))

        if all(isinstance(v, (float, int)) or v is None for v in values):
            values = [None if v is None else v for v in values]
            if all(isinstance(v, int) for v in values if not np.isnan(v)):
                if any(np.isnan(v) for v in values):
                    return np.array(values, dtype=np.float32)
                return np.array(values, dtype=np.int32)
            return np.array(values, dtype=np.float32)

        if all(isinstance(v, str) or v is None or not v for v in values):
            values = [None if v is None or not v else v for v in values]
            return np.array(values, dtype=h5py.string_dtype(encoding="utf-8"))

        try:
            values = [json.dumps(v) for v in values]
            return np.array(values, dtype=h5py.string_dtype(encoding="utf-8"))
        except Exception:
            raise ValueError(
                "Unsupported data types in 'values'. Expected all int, all float, or all str."
            )

    def save(self, filename):
        """Save the data to an HDF5 file."""
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
