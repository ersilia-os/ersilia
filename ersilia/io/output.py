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
from ..default import FEATURE_MERGE_PATTERN
from ..db.hubdata.interfaces import AirtableInterface


class DataFrame(object):
    def __init__(self, data, columns):
        self.data = data
        self.columns = columns

    @staticmethod
    def _is_h5(file_name):
        extension = file_name.split(".")[-1]
        if extension == "h5":
            return True
        else:
            return False

    @staticmethod
    def _get_delimiter(file_name):
        extension = file_name.split(".")[-1]
        if extension == "tsv":
            return "\t"
        else:
            return ","

    def decompose(self):
        features = self.columns[2:]
        keys = [r[0] for r in self.data]
        inputs = [r[1] for r in self.data]
        values = [r[2:] for r in self.data]
        return {"keys": keys, "inputs": inputs, "features": features, "values": values}

    def write_hdf5(self, file_name):
        res = self.decompose()
        hdf5 = Hdf5Data(
            values=res["values"],
            keys=res["keys"],
            inputs=res["inputs"],
            features=res["features"],
        )
        hdf5.save(file_name)

    def write_text(self, file_name, delimiter=None):
        with open(file_name, "w", newline="") as f:
            if delimiter is None:
                delimiter = self._get_delimiter(file_name)
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(self.columns)
            for i, row in enumerate(self.data):
                writer.writerow(row)

    def write(self, file_name, delimiter=None):
        if self._is_h5(file_name):
            self.write_hdf5(file_name)
        else:
            self.write_text(file_name, delimiter=delimiter)


class ResponseRefactor(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        self._expect_meta = None
        self._meta = None

    def _has_meta(self, result):
        if self._expect_meta is not None:
            return self._expect_meta
        try:
            r = result["result"]
            m = result["meta"]
            self._expect_meta = True
        except:
            self._expect_meta = False

    def _get_result(self, result):
        if self._expect_meta is None:
            self._has_meta(result)
        if self._expect_meta:
            return result["result"]
        else:
            return result

    def _get_meta(self, result):
        if self._meta is not None:
            return self._meta
        if self._expect_meta is None:
            self._has_meta(result)
        if self._expect_meta:
            return result["meta"]
        else:
            return None

    def _nullify_meta(self, meta, result):
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

    def refactor_response(self, result):
        r = self._get_result(result)
        m = self._get_meta(result)
        m = self._nullify_meta(m, r)
        self._meta = m
        return r


class GenericOutputAdapter(ResponseRefactor):
    def __init__(self, config_json):
        ResponseRefactor.__init__(self, config_json=config_json)
        self.api_schema = None
        self._schema = None
        self._array_types = set(
            ["array", "numeric_array", "string_array", "mixed_array"]
        )

    @staticmethod
    def _is_string(output):
        if type(output) is str:
            return True
        else:
            return False

    @staticmethod
    def _extension(filename):
        return filename.split(".")[-1]

    def _has_extension(self, output, extension):
        if not self._is_string(output):
            return False
        ext = output.split(".")[-1]
        if ext == extension:
            return True
        else:
            return False

    def __pure_dtype(self, k):
        t = self._schema[k]["type"]
        return t

    def __array_shape(self, k):
        s = self._schema[k]["shape"]
        return s[0]  # TODO work with tensors

    def __meta_by_key(self, k):
        return self._schema[k]["meta"]

    def __cast_values(self, vals, dtypes, output_keys):
        v = []
        for v_, t_, k_ in zip(vals, dtypes, output_keys):
            if t_ in self._array_types:
                if v_ is None:
                    v_ = [None] * self.__array_shape(k_)
                v += v_
            else:
                v += [v_]
        return v

    def _guess_pure_dtype_if_absent(self, vals):
        pdt = PureDataTyper(vals)
        dtype = pdt.get_type()
        self.logger.debug("Guessed pure datatype: {0}".format(dtype))
        return dtype["type"]

    def __expand_output_keys(self, vals, output_keys):
        output_keys_expanded = []
        if len(output_keys) == 1:
            merge_key = False
        else:
            merge_key = True
        current_pure_dtype = {}
        for v, ok in zip(vals, output_keys):
            self.logger.debug("Data: {0}".format(ok))
            self.logger.debug("Values: {0}".format(v))
            m = self.__meta_by_key(ok)
            if ok not in current_pure_dtype:
                t = self.__pure_dtype(ok)
                if t is None:
                    t = self._guess_pure_dtype_if_absent(v)
                current_pure_dtype[ok] = t
            else:
                t = current_pure_dtype[ok]
            self.logger.debug("Pure datatype: {0}".format(t))
            if t in self._array_types:
                assert m is not None
                if v is not None:
                    assert len(m) == len(v)
                if merge_key:
                    output_keys_expanded += [
                        "{0}{1}{2}".format(ok, FEATURE_MERGE_PATTERN, m_) for m_ in m
                    ]
                else:
                    output_keys_expanded += ["{0}".format(m_) for m_ in m]
            else:
                output_keys_expanded += [ok]
        return output_keys_expanded

    def _get_outputshape_from_airtable(self, model_id):
        airtable_interface = AirtableInterface(config_json=self.config_json)
        output_shape = " "
        for record in airtable_interface.items():
            model_idi = record["fields"]["Identifier"]
            try:
                if model_idi == model_id:
                    output_shape = record["fields"]["Output Shape"]
            except KeyError:
                self.logger.warning("The Output Shape field is empty")
        return output_shape

    def _to_dataframe(self, result, model_id):
        output_shape = self._get_outputshape_from_airtable(model_id)
        result = json.loads(result)
        R = []
        output_keys = None
        output_keys_expanded = None
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
                dtypes = [self.__pure_dtype(k) for k in output_keys]
                if output_keys_expanded is None:
                    output_keys_expanded = self.__expand_output_keys(vals, output_keys)
                vals = self.__cast_values(vals, dtypes, output_keys)
            R += [[inp["key"], inp["input"]] + vals]
        columns = ["key", "input"] + output_keys_expanded
        df = DataFrame(data=R, columns=columns)
        return df

    def meta(self):
        if self._meta is None:
            self.logger.error(
                "Meta not available, run some adapations first and it will be inferred atomatically"
            )
        else:
            return self._meta

    def merge(self, subfiles, output_file):
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

    def adapt(self, result, output, model_id=None, api_name=None):
        if model_id is not None and api_name is not None and self.api_schema is None:
            self.api_schema = ApiSchema(model_id=model_id, config_json=self.config_json)
        if self.api_schema is not None:
            if self.api_schema.isfile():
                self._schema = self.api_schema.get_output_by_api(api_name)
        else:
            self.api_schema = None
        if output is not None and self._schema is None:
            raise Exception
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


class DictlistDataframeConverter(GenericOutputAdapter):
    def __init__(self, config_json):
        GenericOutputAdapter.__init__(self, config_json=config_json)

    def dictlist2dataframe(self, dl, model_id, api_name):
        tmp_dir = tempfile.mkdtemp(prefix="ersilia-")
        df_file = os.path.join(tmp_dir, "data.csv")
        self.adapt(dl, df_file, model_id, api_name)
        df = Dataframe()
        df.from_csv(df_file)
        return df

    def __nan_to_none(self, x):
        if np.isnan(x):
            return None
        return x

    def dataframe2dictlist(self, df, model_id, api_name):
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
    def __init__(self, file_names):
        ft = FileTyper(file_names[0])
        self.is_hdf5 = ft.is_hdf5()
        self.file_names = file_names

    def stack_text(self, output):
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

    def stack_hdf5(self, output):
        stacker = Hdf5DataStacker(self.file_names)
        stacker.stack(output)

    def stack(self, output):
        if self.is_hdf5:
            self.stack_hdf5(output)
        else:
            self.stack_text(output)
