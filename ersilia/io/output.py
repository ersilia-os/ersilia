import collections
import csv
import json
import os
import random
import time
from itertools import islice

import numpy as np

from .. import ErsiliaBase
from ..default import (
    FEATURE_MERGE_PATTERN,
)
from ..serve.schema import ApiSchema
from ..utils.exceptions_utils.api_exceptions import UnprocessableInputError
from ..utils.logging import logger, make_temp_dir
from .dataframe import Dataframe
from .readers.file import FileTyper


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
        self.logger = logger
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

    def write_hdf5(self, file_name: str, append: bool = False):
        """
        Writes the DataFrame to an HDF5 file.

        Parameters
        ----------
        file_name : str
            The name of the file to write to.
        append : bool, optional
            If True, append rows to an existing file written by this method.
        """
        from ..utils.hdf5 import Hdf5Data

        res = self.decompose()
        hdf5 = Hdf5Data(
            values=res["values"],
            keys=res["keys"],
            inputs=res["inputs"],
            features=res["features"],
            dtype=self.dtype,
            dim=self.dim,
        )

        hdf5.save(file_name, append=append)

    def write_text(
        self,
        file_name: str,
        delimiter: str = None,
        chunksize: int = 50_000,
        append: bool = False,
    ):
        """
        Writes the DataFrame to a text file, wrapping any string-valued field in quotes,
        except for the first and second columns, which are never quoted.
        If append is True, rows are appended to an existing file without a header.
        """
        if delimiter is None:
            delimiter = self._get_delimiter(file_name)

        cols = list(self.columns)
        ncols = len(cols)

        it = iter(self.data)

        t0 = time.perf_counter()
        mode = "a" if append else "w"
        with open(file_name, mode, newline="", buffering=1024 * 1024) as f:
            writer = csv.writer(
                f,
                delimiter=delimiter,
                lineterminator="\n",
                quoting=csv.QUOTE_MINIMAL,
                quotechar='"',
                doublequote=True,
            )
            if not append:
                writer.writerow(cols)
            self.logger.debug(
                f"write_text header dt={(time.perf_counter() - t0):.6f}s file={file_name}"
            )

            wrote = 0
            while True:
                t = time.perf_counter()
                chunk = list(islice(it, chunksize))
                self.logger.debug(
                    f"write_text chunk_take rows={len(chunk)} dt={(time.perf_counter() - t):.6f}s"
                )
                if not chunk:
                    break

                # csv.writer already writes None as an empty field and other
                # values as str(v), so rows only need fixing when the length is off.
                t = time.perf_counter()
                out_rows = [
                    row
                    if len(row) == ncols
                    else (list(row) + [None] * (ncols - len(row)))[:ncols]
                    for row in chunk
                ]
                self.logger.debug(
                    f"write_text chunk_format rows={len(out_rows)} dt={(time.perf_counter() - t):.6f}s"
                )

                t = time.perf_counter()
                writer.writerows(out_rows)
                wrote += len(out_rows)
                self.logger.debug(
                    f"write_text chunk_write rows={len(out_rows)} dt={(time.perf_counter() - t):.6f}s"
                )

        self.logger.debug(
            f"write_text done rows={wrote} dt_total={(time.perf_counter() - t0):.6f}s file={file_name}"
        )

    def write(
        self,
        file_name: str,
        delimiter: str = None,
        append: bool = False,
        check_unprocessable: bool = True,
    ):
        """
        Writes the DataFrame to a file, determining the format based on the file extension.

        Parameters
        ----------
        file_name : str
            The name of the file to write to.
        delimiter : str, optional
            The delimiter to use in the text file (default is None).
        append : bool, optional
            If True, append rows to a file previously written by this method.
        check_unprocessable : bool, optional
            If True, raise when the data is a single unprocessable input. Disable
            when writing one chunk of a larger output.
        """
        if check_unprocessable and self._is_unprocessable_input():
            raise UnprocessableInputError()
        if self._is_h5(file_name):
            self.write_hdf5(file_name, append=append)
        else:
            self.write_text(file_name, delimiter=delimiter, append=append)


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

    def _to_dataframe(self, result) -> DataFrame:
        t0 = time.perf_counter()
        result = self._try_serialize_str_to_json(result)
        t1 = time.perf_counter()
        logger.debug(f"_to_dataframe serialize dt={(t1 - t0):.6f}s")

        in_cols = self.input_columns
        out_cols = self.output_columns
        cols = in_cols + out_cols
        dim = len(cols)
        cast = self.output_dtype

        data = []
        append_row = data.append

        t2 = time.perf_counter()
        for r in result:
            rin = r["input"]
            rout = r["output"]

            row = [rin[c] for c in in_cols]
            for c in out_cols:
                v = rout[c]
                row.append(cast(v) if v is not None else None)

            append_row(row)

        t3 = time.perf_counter()
        logger.debug(f"_to_dataframe build_rows rows={len(data)} dt={(t3 - t2):.6f}s")

        df = DataFrame(data=data, columns=cols, dtype=cast, dim=dim)
        t4 = time.perf_counter()
        logger.debug(
            f"GenericOutputAdapter_to_dataframe construct_df dt={(t4 - t3):.6f}s dt_total={(t4 - t0):.6f}s"
        )
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
        return result, df

    def write_chunk(self, result: list, output: str, append: bool):
        """
        Writes one chunk of standardized results to a CSV, TSV or HDF5 file.

        Used to stream large runs to disk batch by batch instead of holding
        the full output in memory.

        Parameters
        ----------
        result : list
            Standardized results ({"input": ..., "output": ...} dicts).
        output : str
            The output file name.
        append : bool
            If True, append to the file written by a previous chunk.
        """
        df = self._to_dataframe(result)
        delimiter = "\t" if self._has_extension(output, "tsv") else ","
        df.write(output, delimiter=delimiter, append=append, check_unprocessable=False)

    @staticmethod
    def _row_values(result, n):
        # Same positional mapping as StandardCSVRunApi._standardize_output: the
        # first value is the input, the next n are the output columns in order.
        if result is None:
            return [None] * n
        values = result if isinstance(result, list) else list(result.values())
        values = values[1 : n + 1]
        if len(values) < n:
            values += [None] * (n - len(values))
        return values

    @staticmethod
    def _float_rows(values):
        # Convert the whole batch at once; rows holding None are redone cell by
        # cell because numpy would turn None into nan instead of an empty field.
        try:
            rows = np.array(values, dtype=np.float64).tolist()
        except (TypeError, ValueError):
            return [[None if v is None else float(v) for v in row] for row in values]
        for i, row in enumerate(values):
            if None in row:
                rows[i] = [None if v is None else float(v) for v in row]
        return rows

    def write_batch(self, inputs: list, results: list, output: str, append: bool):
        """
        Writes one batch of merged API results to a CSV, TSV or HDF5 file.

        Produces the same file as standardizing the batch and calling
        write_chunk, but builds the rows in one pass and converts float
        outputs as a single matrix instead of cell by cell.

        Parameters
        ----------
        inputs : list
            Input records ({"key": ..., "input": ...}), one per result.
        results : list
            One dict per input: the input first, then the output values in
            column order.
        output : str
            The output file name.
        append : bool
            If True, append to the file written by a previous batch.
        """
        n = len(self.output_columns)
        cast = self.output_dtype
        cols = self.input_columns + self.output_columns
        keys = [d["key"] for d in inputs]
        smiles = [d["input"] for d in inputs]
        values = [self._row_values(r, n) for r in results]

        if self._has_extension(output, "h5"):
            from ..utils.hdf5 import Hdf5Data

            if cast is float:
                try:
                    # None becomes nan, the fill value Hdf5Data uses for floats.
                    values = np.array(values, dtype=np.float64)
                except (TypeError, ValueError):
                    pass
            hdf5 = Hdf5Data(
                values=values,
                keys=keys,
                inputs=smiles,
                features=self.output_columns,
                dtype=cast,
                dim=len(cols),
            )
            hdf5.save(output, append=append)
            return

        if cast is float:
            rows = self._float_rows(values)
        else:
            rows = [[None if v is None else cast(v) for v in row] for row in values]
        data = [[k, s, *row] for k, s, row in zip(keys, smiles, rows)]
        df = DataFrame(data=data, columns=cols, dtype=cast, dim=len(cols))
        delimiter = "\t" if self._has_extension(output, "tsv") else ","
        df.write_text(output, delimiter=delimiter, append=append)

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
        adapted_result, _ = self._adapt_generic(result, output, model_id, api_name)
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
        from ..utils.hdf5 import Hdf5DataStacker

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
