import collections
import json
import os
import random

import numpy as np

from .. import ErsiliaBase
from ..default import (
    FEATURE_MERGE_PATTERN,
)
from ..serve.schema import ApiSchema
from ..utils.exceptions_utils.api_exceptions import UnprocessableInputError
from ..utils.hdf5 import Hdf5Data, Hdf5DataStacker
from ..utils.logging import make_temp_dir
from .dataframe import Dataframe
from .readers.file import FileTyper


class DataFrame(object):
    """
    A class used to represent a DataFrame.

    Attributes
    ----------
    data : list
        The data of the DataFrame.
    columns : list
        The column names of the DataFrame.
    dtype : any
        Data type for the output part of the data
    dim: int
        The number of dimensions of the data.

    Methods
    -------
    decompose()
        Decomposes the DataFrame into its components.
    write_hdf5(file_name)
        Writes the DataFrame to an HDF5 file.
    write_text(file_name, delimiter=None)
        Writes the DataFrame to a text file.
    write(file_name, delimiter=None)
        Writes the DataFrame to a file, determining the format based on the file extension.
    """

    def __init__(self, data: list, columns: list, dtype: any, dim: int):
        self.data = data
        self.dim = dim
        self.dtype = dtype
        self.columns = columns

    def _is_unprocessable_input(self) -> bool:
        if len(self.data) != 1:
            return False
        row = self.data[0]
        return row[0] == "UNPROCESSABLE_INPUT" and row[1] == "UNPROCESSABLE_INPUT"

    @staticmethod
    def _is_h5(file_name: str) -> bool:
        """
        Checks if the file is an HDF5 file.

        Parameters
        ----------
        file_name : str
            The name of the file.

        Returns
        -------
        bool
            True if the file is an HDF5 file, False otherwise.
        """
        extension = file_name.split(".")[-1]
        if extension == "h5":
            return True
        else:
            return False

    @staticmethod
    def _get_delimiter(file_name: str) -> str:
        """
        Gets the delimiter based on the file extension.

        Parameters
        ----------
        file_name : str
            The name of the file.

        Returns
        -------
        str
            The delimiter.
        """
        extension = file_name.split(".")[-1]
        if extension == "tsv":
            return "\t"
        else:
            return ","

    def decompose(self) -> dict:
        """
        Decomposes the DataFrame into its components.

        Returns
        -------
        dict
            A dictionary containing keys, inputs, features, and values.
        """
        features = self.columns[2:]
        keys = [r[0] for r in self.data]
        inputs = [r[1] for r in self.data]
        values = [r[2:] for r in self.data]
        return {"keys": keys, "inputs": inputs, "features": features, "values": values}

    def write_hdf5(self, file_name: str):
        """
        Writes the DataFrame to an HDF5 file.

        Parameters
        ----------
        file_name : str
            The name of the file to write to.
        """
        res = self.decompose()
        hdf5 = Hdf5Data(
            values=res["values"],
            keys=res["keys"],
            inputs=res["inputs"],
            features=res["features"],
            dtype=self.dtype,
            dim=self.dim,
        )

        hdf5.save(file_name)

    def write_text(self, file_name: str, delimiter: str = None):
        """
        Writes the DataFrame to a text file, wrapping any string-valued field in quotes,
        except for the first and second columns, which are never quoted.
        """
        if delimiter is None:
            delimiter = self._get_delimiter(file_name)

        none_str = ""
        with open(file_name, "w", newline="") as f:
            f.write(delimiter.join(self.columns) + "\n")

            for row in self.data:
                out_fields = []
                for idx, val in enumerate(row):
                    if val is None:
                        text = none_str
                    else:
                        text = str(val)
                        if isinstance(val, str) and idx > 1:
                            text = f'"{text}"'
                    out_fields.append(text)

                f.write(delimiter.join(out_fields) + "\n")

    def write(self, file_name: str, delimiter: str = None):
        """
        Writes the DataFrame to a file, determining the format based on the file extension.

        Parameters
        ----------
        file_name : str
            The name of the file to write to.
        delimiter : str, optional
            The delimiter to use in the text file (default is None).
        """
        if self._is_unprocessable_input():
            raise UnprocessableInputError()
        if self._is_h5(file_name):
            self.write_hdf5(file_name)
        else:
            self.write_text(file_name, delimiter=delimiter)


class ResponseRefactor(ErsiliaBase):
    """
    A class used to refactor API responses.

    Attributes
    ----------
    config_json : dict
        Configuration settings in JSON format.

    Methods
    -------
    refactor_response(result)
        Refactors the API response.
    """

    def __init__(self, config_json: dict):
        """
        Parameters
        ----------
        config_json : dict
            Configuration settings in JSON format.
        """
        ErsiliaBase.__init__(self, config_json=config_json)
        self._expect_meta = None
        self._meta = None

    def _has_meta(self, result: dict) -> bool:
        """
        Checks if the result has meta information.

        Parameters
        ----------
        result : dict
            The result to check.

        Returns
        -------
        bool
            True if the result has meta information, False otherwise.
        """
        if self._expect_meta is not None:
            return self._expect_meta
        try:
            r = result["result"]  # noqa: F841
            m = result["meta"]  # noqa: F841
            self._expect_meta = True
        except:
            self._expect_meta = False

    def _get_result(self, result: dict) -> dict:
        """
        Gets the result from the response.

        Parameters
        ----------
        result : dict
            The response to get the result from.

        Returns
        -------
        dict
            The result.
        """
        if self._expect_meta is None:
            self._has_meta(result)
        if self._expect_meta:
            return result["result"]
        else:
            return result

    def _get_meta(self, result: dict) -> dict:
        """
        Gets the meta information from the response.

        Parameters
        ----------
        result : dict
            The response to get the meta information from.

        Returns
        -------
        dict
            The meta information.
        """
        if self._meta is not None:
            return self._meta
        if self._expect_meta is None:
            self._has_meta(result)
        if self._expect_meta:
            return result["meta"]
        else:
            return None

    def _nullify_meta(self, meta: dict, result: dict) -> dict:
        """
        Nullifies the meta information.

        Parameters
        ----------
        meta : dict
            The meta information to nullify.
        result : dict
            The result to use for nullifying the meta information.

        Returns
        -------
        dict
            The nullified meta information.
        """
        m = {}
        one_output = random.choice(result)
        for k, v in one_output.items():
            if meta is None:
                m[k] = None
            else:
                if k not in meta:
                    m[k] = None
                else:
                    m[k] = meta[k]
        return m

    def refactor_response(self, result: dict) -> dict:
        """
        Refactors the API response.

        Parameters
        ----------
        result : dict
            The API response to refactor.

        Returns
        -------
        dict
            The refactored API response.
        """
        r = self._get_result(result)
        m = self._get_meta(result)
        m = self._nullify_meta(m, r)
        self._meta = m
        return r


class GenericOutputAdapter(ResponseRefactor):
    """
    A class used to adapt generic outputs.

    Attributes
    ----------
    model_id : str
        The model identifier.
    config_json : dict
        Configuration settings in JSON format.

    Methods
    -------
    adapt(result, output, model_id=None, api_name=None)
        Adapts the output based on the result and model.
    """

    def __init__(
        self, model_id: str = None, columns_info: dict = None, config_json: dict = None
    ):
        ResponseRefactor.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.columns_info = columns_info
        self.input_columns = ["key", "input"]
        self.output_columns = self.columns_info["name"]
        self.output_dtype = self._get_output_dtype()

    def _get_output_dtype(self):
        dtypes = sorted(set(self.columns_info["type"]))
        if dtypes == ["string"]:
            return str
        elif dtypes == ["float"]:
            return float
        elif dtypes == ["integer"]:
            return int
        else:
            return float

    def _try_serialize_str_to_json(self, data):
        if self._is_string(data):
            return json.loads(data)
        return data

    @staticmethod
    def _is_string(output: any) -> bool:
        if type(output) is str:
            return True
        else:
            return False

    @staticmethod
    def _extension(filename: str) -> str:
        return filename.split(".")[-1]

    def _has_extension(self, output: any, extension: str) -> bool:
        if not self._is_string(output):
            return False
        ext = output.split(".")[-1]
        if ext == extension:
            return True
        else:
            return False

    def _to_dataframe(self, result: dict) -> DataFrame:
        result = self._try_serialize_str_to_json(result)
        input_data = collections.defaultdict(list)
        output_data = collections.defaultdict(list)
        for r in result:
            for c in self.input_columns:
                input_data[c] += [r["input"][c]]
            for c in self.output_columns:
                output_data[c] += [self.output_dtype(r["output"][c])]
        data = []
        columns = self.input_columns + self.output_columns
        dtype = self.output_dtype

        R = []
        for c in self.input_columns:
            R += [input_data[c]]
        for c in self.output_columns:
            R += [output_data[c]]
        data = list(map(list, zip(*R)))
        df = DataFrame(data=data, columns=columns, dtype=dtype, dim=len(columns))
        return df

    def merge(self, subfiles: list, output_file: str):
        """
        Merges the subfiles into a single output file.

        Parameters
        ----------
        subfiles : list
            The list of subfiles to merge.
        output_file : str
            The name of the output file.
        """
        self.logger.debug(
            "Merging {0} files into {1}".format(len(subfiles), output_file)
        )
        extensions = set([self._extension(x) for x in subfiles + [output_file]])
        assert len(extensions) == 1
        if self._has_extension(output_file, "json"):
            data = []
            for subfile in subfiles:
                with open(subfile, "r") as f:
                    data += json.load(f)
            with open(output_file, "w") as f:
                json.dump(data, f, indent=4)
        else:
            with open(output_file, "w") as fo:
                use_header = True
                for subfile in subfiles:
                    with open(subfile, "r") as fi:
                        if not use_header:
                            next(fi)
                        for l in fi:
                            fo.write(l)
                    use_header = False

    def _adapt_generic(
        self, result: str, output: str, model_id: str = None, api_name: str = None
    ) -> dict:
        if api_name != "run":
            self.logger.debug("Api was not run")
            raise Exception(
                "Only 'run' api is supported for now... Please open an issue if you feel Ersilia needs to add"
            )
        if model_id is None:
            self.logger.debug("Model ID is None")
            pass
        if output is None:
            self.logger.debug("Output is None")
            return None
        if self._has_extension(output, "csv"):
            extension = "csv"
        elif self._has_extension(output, "tsv"):
            extension = "tsv"
        elif self._has_extension(output, "h5"):
            extension = "h5"
        elif self._has_extension(output, "json"):
            extension = "json"
        else:
            extension = None
        self.logger.debug(f"Extension: {extension}")
        df = self._to_dataframe(result)
        delimiters = {"csv": ",", "tsv": "\t"}
        if extension in ["tsv", "csv"]:
            df.write(output, delimiter=delimiters[extension])
        elif extension in ["h5"]:
            df.write(output)
        elif extension == "json":
            data = json.loads(result)
            with open(output, "w") as f:
                json.dump(data, f, indent=4)
        else:
            pass
        return result

    def adapt(
        self, result: str, output: str, model_id: str = None, api_name: str = None
    ) -> dict:
        """
        Adapts the output based on the result and model.

        Parameters
        ----------
        result : str
            The result to adapt, stringyfied JSON.
        output : str
            The output file name.
        model_id : str, optional
            The model identifier (default is None).
        api_name : str, optional
            The API name (default is None).

        Returns
        -------
        dict
            The adapted result.
        """
        adapted_result = self._adapt_generic(result, output, model_id, api_name)
        return adapted_result


class DictlistDataframeConverter(GenericOutputAdapter):
    """
    A class used to convert between dictlist and DataFrame.

    Attributes
    ----------
    config_json : dict
        Configuration settings in JSON format.

    Methods
    -------
    dictlist2dataframe(dl, model_id, api_name)
        Converts a dictlist to a DataFrame.
    dataframe2dictlist(df, model_id, api_name)
        Converts a DataFrame to a dictlist.
    """

    def __init__(self, config_json: dict):
        """
        Parameters
        ----------
        config_json : dict
            Configuration settings in JSON format.
        """
        GenericOutputAdapter.__init__(self, config_json=config_json)

    def dictlist2dataframe(self, dl: dict, model_id: str, api_name: str) -> Dataframe:
        """
        Converts a dictlist to a DataFrame.

        Parameters
        ----------
        dl : dict
            The dictlist to convert.
        model_id : str
            The model identifier.
        api_name : str
            The API name.

        Returns
        -------
        Dataframe
            The converted DataFrame.
        """
        tmp_dir = make_temp_dir(prefix="ersilia-")
        df_file = os.path.join(tmp_dir, "data.csv")
        self.adapt(dl, df_file, model_id, api_name)
        df = Dataframe()
        df.from_csv(df_file)
        return df

    def __nan_to_none(self, x: any) -> any:
        """
        Converts NaN values to None.

        Parameters
        ----------
        x : any
            The value to convert.

        Returns
        -------
        any
            The converted value.
        """
        if np.isnan(x):
            return None
        return x

    def dataframe2dictlist(self, df: DataFrame, model_id: str, api_name: str) -> list:
        """
        Converts a DataFrame to a dictlist.

        Parameters
        ----------
        df : DataFrame
            The DataFrame to convert.
        model_id : str
            The model identifier.
        api_name : str
            The API name.

        Returns
        -------
        list
            The converted dictlist.
        """
        schema = ApiSchema(
            model_id=model_id, config_json=self.config_json
        ).get_output_by_api(api_name=api_name)
        the_keys = [k for k, _ in schema.items()]
        if len(the_keys) == 1:
            the_key = the_keys[0]
        else:
            the_key = None
        result = []
        features = df.features
        grouped_features_idxs = collections.defaultdict(list)
        grouped_features = collections.defaultdict(list)
        for i, f in enumerate(features):
            if the_key is None:
                splitted = f.split(FEATURE_MERGE_PATTERN)
                if len(splitted) == 2:
                    g, f = f.split(FEATURE_MERGE_PATTERN)
                else:
                    g = f
            else:
                g = the_key
            grouped_features_idxs[g] += [i]
            grouped_features[g] += [f]
        # Reorder to match schema, just to be sure
        for k, v in grouped_features.items():
            ords = dict((k_, i_) for i_, k_ in enumerate(v))
            if schema[k]["meta"] is not None:
                ord_idxs = [ords[v_] for v_ in schema[k]["meta"]]
                grouped_features[k] = [v[idx] for idx in ord_idxs]
                w = grouped_features_idxs[k]
                grouped_features_idxs[k] = [w[idx] for idx in ord_idxs]
        for r in df.iterrows():
            output = {}
            for k, idxs in grouped_features_idxs.items():
                v = [self.__nan_to_none(x) for x in r["values"][idxs].tolist()]
                if len(v) == 1:
                    v = v[0]
                output[k] = v
            res = {
                "input": {"key": r["key"], "input": r["input"], "text": None},
                "output": output,
            }
            result += [res]
        return result


class TabularOutputStacker(object):
    """
    A class used to stack tabular outputs.

    Attributes
    ----------
    file_names : list
        The list of file names to stack.

    Methods
    -------
    stack(output)
        Stacks the files into a single output file.
    """

    def __init__(self, file_names: list):
        """
        Parameters
        ----------
        file_names : list
            The list of file names to stack.
        """
        ft = FileTyper(file_names[0])
        self.is_hdf5 = ft.is_hdf5()
        self.file_names = file_names

    def stack_text(self, output: str):
        """
        Stacks text files into a single output file.

        Parameters
        ----------
        output : str
            The name of the output file.
        """
        has_header = False
        with open(output, "w") as f0:
            for fn in self.file_names:
                with open(fn, "r") as f1:
                    header = next(f1)
                    if not has_header:
                        f0.write(header)
                        has_header = True
                    for l in f1:
                        f0.write(l)

    def stack_hdf5(self, output: str):
        """
        Stacks HDF5 files into a single output file.

        Parameters
        ----------
        output : str
            The name of the output file.
        """
        stacker = Hdf5DataStacker(self.file_names)
        stacker.stack(output)

    def stack(self, output: str):
        """
        Stacks the files into a single output file.

        Parameters
        ----------
        output : str
            The name of the output file.
        """
        if self.is_hdf5:
            self.stack_hdf5(output)
        else:
            self.stack_text(output)
