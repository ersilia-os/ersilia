
import csv
import os
import json
import random
import numpy as np
import tempfile
import collections
from .dataframe import Dataframe
from .readers.file import FileTyper
from .pure import PureDataTyper
from ..serve.schema import ApiSchema
from .. import ErsiliaBase
from ..utils.hdf5 import Hdf5Data, Hdf5DataStacker
from ..db.hubdata.interfaces import JsonModelsInterface
from ..default import FEATURE_MERGE_PATTERN, PACK_METHOD_FASTAPI
from ..utils.paths import resolve_pack_method
from ..utils.logging import make_temp_dir


class DataFrame(object):
    """
    A class used to represent a DataFrame.

    Attributes
    ----------
    data : list
        The data of the DataFrame.
    columns : list
        The column names of the DataFrame.

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

    def __init__(self, data: list, columns: list):
        """
        Parameters
        ----------
        data : list
            The data of the DataFrame.
        columns : list
            The column names of the DataFrame.
        """
        self.data = data
        self.columns = columns

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
        )
        hdf5.save(file_name)

    def write_text(self, file_name: str, delimiter: str = None):
        """
        Writes the DataFrame to a text file.

        Parameters
        ----------
        file_name : str
            The name of the file to write to.
        delimiter : str, optional
            The delimiter to use in the text file (default is None).
        """
        with open(file_name, "w", newline="") as f:
            if delimiter is None:
                delimiter = self._get_delimiter(file_name)
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(self.columns)
            for i, row in enumerate(self.data):
                writer.writerow(row)

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
            r = result["result"]
            m = result["meta"]
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

    def __init__(self, model_id: str = None, config_json: dict = None):
        """
        Parameters
        ----------
        model_id : str, optional
            The model identifier (default is None).
        config_json : dict, optional
            Configuration settings in JSON format (default is None).
        """
        ResponseRefactor.__init__(self, config_json=config_json)
        self.api_schema = None
        self._schema = None
        self._array_types = set(
            ["array", "numeric_array", "string_array", "mixed_array"]
        )
        self.model_id = model_id
        self.was_fast_api = (
            resolve_pack_method(self._get_bundle_location(self.model_id))
            == PACK_METHOD_FASTAPI
        )

    @staticmethod
    def _is_string(output: any) -> bool:
        """
        Checks if the output is a string.

        Parameters
        ----------
        output : any
            The output to check.

        Returns
        -------
        bool
            True if the output is a string, False otherwise.
        """
        if type(output) is str:
            return True
        else:
            return False

    @staticmethod
    def _extension(filename: str) -> str:
        """
        Gets the file extension.

        Parameters
        ----------
        filename : str
            The name of the file.

        Returns
        -------
        str
            The file extension.
        """
        return filename.split(".")[-1]

    def _has_extension(self, output: any, extension: str) -> bool:
        """
        Checks if the output has the specified extension.

        Parameters
        ----------
        output : any
            The output to check.
        extension : str
            The extension to check for.

        Returns
        -------
        bool
            True if the output has the specified extension, False otherwise.
        """
        if not self._is_string(output):
            return False
        ext = output.split(".")[-1]
        if ext == extension:
            return True
        else:
            return False

    def __pure_dtype(self, k: str) -> str:
        """
        Gets the pure data type for the specified key.

        Parameters
        ----------
        k : str
            The key to get the pure data type for.

        Returns
        -------
        str
            The pure data type.
        """
        t = self._schema[k]["type"]
        return t

    def __array_shape(self, k: str) -> int:
        """
        Gets the array shape for the specified key.

        Parameters
        ----------
        k : str
            The key to get the array shape for.

        Returns
        -------
        int
            The array shape.
        """
        s = self._schema[k]["shape"]
        return s[0]  # TODO work with tensors

    def __meta_by_key(self, k: str) -> dict:
        """
        Gets the meta information for the specified key.

        Parameters
        ----------
        k : str
            The key to get the meta information for.

        Returns
        -------
        dict
            The meta information.
        """
        return self._schema[k]["meta"]

    def __cast_values(self, vals: list, dtypes: list, output_keys: list) -> list:
        """
        Casts the values to the specified data types.

        Parameters
        ----------
        vals : list
            The values to cast.
        dtypes : list
            The data types to cast to.
        output_keys : list
            The output keys.

        Returns
        -------
        list
            The casted values.
        """
        v = []
        for v_, t_, k_ in zip(vals, dtypes, output_keys):
            self.logger.debug(v_)
            self.logger.debug(t_)
            self.logger.debug(k_)
            if t_ in self._array_types:
                if v_ is None:
                    v_ = [None] * self.__array_shape(k_)
                v += v_
            else:
                v += [v_]
        return v

    def _guess_pure_dtype_if_absent(self, vals: list) -> str:
        """
        Guesses the pure data type if it is absent.

        Parameters
        ----------
        vals : list
            The values to guess the data type for.

        Returns
        -------
        str
            The guessed pure data type.
        """
        pdt = PureDataTyper(vals, model_id=self.model_id)
        dtype = pdt.get_type()
        self.logger.debug("Guessed pure datatype: {0}".format(dtype))
        if dtype is None:
            return None
        else:
            return dtype["type"]

    def __expand_output_keys(self, vals: list, output_keys: list) -> list:
        """
        Expands the output keys.

        Parameters
        ----------
        vals : list
            The values to expand the output keys for.
        output_keys : list
            The output keys to expand.

        Returns
        -------
        list
            The expanded output keys.
        """
        output_keys_expanded = []
        if len(output_keys) == 1:
            merge_key = False
        else:
            merge_key = True
        current_pure_dtype = {}
        for key_index, val_out in enumerate(zip(vals, output_keys)):
            v, ok = val_out
            self.logger.debug("Data: {0}".format(ok))
            self.logger.debug("Values: {0}".format(str(v)[:100]))
            m = self.__meta_by_key(ok)
            if ok not in current_pure_dtype:
                self.logger.debug("Getting pure dtype for {0}".format(ok))
                if self.dtypes is not None:
                    t = self.dtypes[key_index]
                else:
                    t = self.__pure_dtype(ok)
                self.logger.debug("This is the pure datatype: {0}".format(t))
                if t is None:
                    t = self._guess_pure_dtype_if_absent(v)
                    self.logger.debug("Guessed absent pure datatype: {0}".format(t))
                current_pure_dtype[ok] = t
            else:
                t = current_pure_dtype[ok]
            self.logger.debug("Datatype: {0}".format(t))
            if t in self._array_types:
                self.logger.debug(
                    "Datatype has been matched: {0} over {1}".format(
                        t, self._array_types
                    )
                )
                assert m is not None
                if v is not None:
                    if len(m) > len(v):
                        self.logger.debug(
                            "Metadata {0} is longer than values {1}".format(
                                len(m), len(v)
                            )
                        )
                        v = list(v) + [None] * (len(m) - len(v))
                    assert len(m) == len(v)
                if merge_key:
                    self.logger.debug("Merge key is {0}".format(merge_key))
                    output_keys_expanded += [
                        "{0}{1}{2}".format(ok, FEATURE_MERGE_PATTERN, m_) for m_ in m
                    ]
                else:
                    self.logger.debug("No merge key")
                    output_keys_expanded += ["{0}".format(m_) for m_ in m]
            else:
                output_keys_expanded += [ok]
        return output_keys_expanded

    def _get_outputshape_from_s3_models_json(self, model_id: str) -> str:
        """
        Gets the output shape from the S3 models JSON.

        Parameters
        ----------
        model_id : str
            The model identifier.

        Returns
        -------
        str
            The output shape.
        """
        s3_models_interface = JsonModelsInterface(config_json=self.config_json)
        output_shape = None
        for mdl in s3_models_interface.items():
            _model_id = mdl["Identifier"]
            try:
                if _model_id == model_id:
                    output_shape = mdl["Output Shape"]
            except KeyError:
                self.logger.warning("The Output Shape key is missing")
        return output_shape

    def _get_outputshape(self, model_id: str) -> str:
        """
        Gets the output shape.

        Parameters
        ----------
        model_id : str
            The model identifier.

        Returns
        -------
        str
            The output shape.
        """
        output_shape = self._get_outputshape_from_s3_models_json(model_id)
        if output_shape is None:
            self.logger.warning("Output shape not found")
            output_shape = " "
        return output_shape

    def _to_dataframe(self, result: dict, model_id: str) -> DataFrame:
        """
        Converts the result to a DataFrame.

        Parameters
        ----------
        result : dict
            The result to convert.
        model_id : str
            The model identifier.

        Returns
        -------
        DataFrame
            The converted DataFrame.
        """
        output_shape = self._get_outputshape(model_id)
        result = json.loads(result)
        R = []
        output_keys = None
        output_keys_expanded = None
        self.dtypes = None
        for r in result:
            inp = r["input"]
            out = r["output"]
            if output_shape == "Flexible List":
                vals = [json.dumps(out)]
                output_keys_expanded = ["outcome"]
            else:
                if output_keys is None:
                    output_keys = [k for k in out.keys()]
                vals = [out[k] for k in output_keys]
                # if dtypes has been resolved previously, then it is not necessary to resolve it again
                if self.dtypes is None:
                    self.dtypes = [self.__pure_dtype(k) for k in output_keys]
                are_dtypes_informative = False
                for dtype in self.dtypes:
                    if dtype is not None:
                        are_dtypes_informative = True
                if output_keys_expanded is None:
                    output_keys_expanded = self.__expand_output_keys(vals, output_keys)
                if not are_dtypes_informative:
                    t = self._guess_pure_dtype_if_absent(vals)
                    if len(output_keys) == 1:
                        self.dtypes = [t]
                vals = self.__cast_values(vals, self.dtypes, output_keys)
            R += [[inp["key"], inp["input"]] + vals]
        columns = ["key", "input"] + output_keys_expanded
        df = DataFrame(data=R, columns=columns)
        return df

    def meta(self) -> dict:
        """
        Gets the meta information.

        Returns
        -------
        dict
            The meta information.
        """
        if self._meta is None:
            self.logger.error(
                "Meta not available, run some adapations first and it will be inferred atomatically"
            )
        else:
            return self._meta

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

    def _adapt_generic(self, result: dict, output: str, model_id: str = None, api_name: str = None) -> dict:
        """
        Adapts the output based on the result and model.

        Parameters
        ----------
        result : dict
            The result to adapt.
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
        if model_id is not None and api_name is not None and self.api_schema is None:
            self.api_schema = ApiSchema(model_id=model_id, config_json=self.config_json)
        if self.api_schema is not None:
            if self.api_schema.isfile():
                self._schema = self.api_schema.get_output_by_api(api_name)
        else:
            self.api_schema = None
        if output is not None and self._schema is None:
            raise Exception("Schema not available")
        if self._has_extension(output, "json"):
            data = json.loads(result)
            with open(output, "w") as f:
                json.dump(data, f, indent=4)
        if self._has_extension(output, "csv"):
            df = self._to_dataframe(result, model_id)
            df.write(output)
        if self._has_extension(output, "tsv"):
            df = self._to_dataframe(result, model_id)
            df.write(output, delimiter="\t")
        if self._has_extension(output, "h5"):
            df = self._to_dataframe(result, model_id)
            df.write(output)
        return result

    def _adapt_when_fastapi_was_used(
        self, result: dict, output: str, model_id: str = None, api_name: str = None
    ) -> dict:
        """
        Adapts the output when FastAPI was used.

        Parameters
        ----------
        result : dict
            The result to adapt.
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
        if api_name != "run":
            return None
        if model_id is None:
            return None
        if output is None:
            return None
        if not self.was_fast_api:
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
        delimiters = {"csv": ",", "tsv": "\t"}
        if extension in ["csv", "tsv"]:
            R = []
            for r in json.loads(result):
                inp = r["input"]
                out = r["output"]
                vals = [out[k] for k in out.keys()]
                R += [[inp["key"], inp["input"]] + vals]
            header = ["key", "input"] + [k for k in out.keys()]
            with open(output, "w") as f:
                writer = csv.writer(f, delimiter=delimiters[extension])
                writer.writerow(header)
                for r in R:
                    writer.writerow(r)
        elif extension == "json":
            data = json.loads(result)
            with open(output, "w") as f:
                json.dump(data, f, indent=4)
        elif extension == "h5":
            df = self._to_dataframe(
                result, model_id
            )  # TODO: we can potentially simplify this and get rid of the to_dataframe method for conversion to HDF5.
            df.write(output)
        else:
            pass
        return result

    def adapt(self, result: dict, output: str, model_id: str = None, api_name: str = None) -> dict:
        """
        Adapts the output based on the result and model.

        Parameters
        ----------
        result : dict
            The result to adapt.
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
        adapted_result = self._adapt_when_fastapi_was_used(
            result, output, model_id, api_name
        )
        if adapted_result is None:
            return self._adapt_generic(result, output, model_id, api_name)
        else:
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