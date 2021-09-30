import csv
import json
import random
import itertools
from ..io.pure import PureDataTyper
from .. import ErsiliaBase


class DataFrame(object):
    def __init__(self, data, columns):
        self._is_serialized = False
        self.data = data
        self.columns = columns

    @staticmethod
    def _get_delimiter(file_name):
        extension = file_name.split(".")[-1]
        if extension == "tsv":
            return "\t"
        else:
            return ","

    def serialize(self):
        if self._is_serialized:
            return
        skip_cols = ["key", "input"]
        skip_cols_set = set(skip_cols)
        idx_cols = [
            (j, col) for j, col in enumerate(self.columns) if col not in skip_cols_set
        ]
        new_cols = []
        done_cols = set()
        R = []
        for d in self.data:
            r = []
            for j, col in idx_cols:
                v = [v for v in itertools.chain.from_iterable([d[j]])]
                if col not in done_cols:
                    if len(v) == 1:
                        names = [col]
                    else:
                        names = self._default_feature_names(col, len(v))
                    new_cols += names
                r += v
            R += [r]
        S = []
        for i in range(len(self.data)):
            s = [self.data[i][0], self.data[i][1]] + R[i]
            S += [s]
        self.data = S
        self.columns = skip_cols + new_cols
        self._is_serialized = True

    def write(self, file_name, delimiter=None):
        with open(file_name, "w", newline="") as f:
            if delimiter is None:
                delimiter = self._get_delimiter(file_name)
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(self.columns)
            for i, row in enumerate(self.data):
                writer.writerow(row)



class ResponseRefactor(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.logger.debug("Generic output adapter initialized")
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
        self.logger.debug("Generic output adapter initialized")

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

    def __pure_dtype(self, v):
        dt = PureDataTyper(v)
        t = dt.get_type()
        return t

    def __cast_values(self, v):
        t = self.__pure_dtype(v)
        if t == "array":
            return v
        else:
            return [v]

    @staticmethod
    def __default_feature_names(col, n):
        chars = len(str(n))
        names = []
        for i in range(n):
            i = str(i).zfill(chars)
            names += ["{0}-{1}".format(col, i)]
        return names

    def __expand_output_keys(self, vals, output_keys):
        output_keys_expanded = []
        for v, ok in zip(vals, output_keys):
            m = self._meta[ok]
            t = self.__pure_dtype(v)
            if t == "array":
                if m is None:
                    output_keys_expanded += self._default_feature_names(ok, len(v))
                else:
                    assert len(m) == len(v)
                    output_keys_expanded += ["{0}-{1}".format(ok, m)]
            else:
                output_keys_expanded += [ok]
        return output_keys_expanded

    def _to_dataframe(self, result):
        result = json.loads(result)
        R = []
        output_keys = None
        output_keys_expanded = None
        for r in result:
            inp = r["input"]
            out = r["output"]
            if output_keys is None:
                output_keys = [k for k in out.keys()]
            vals = [out[k] for k in output_keys]
            if output_keys_expanded is None:
                output_keys_expanded = self.__expand_output_keys(vals, output_keys)
            vals = self.__cast_values(vals)
            R += [[inp["key"], inp["input"]] + vals]
        df = DataFrame(data=R, columns=["key", "input"] + output_keys_expanded)
        df.serialize()
        return df

    def meta(self):
        if self._meta is None:
            self.logger.error("Meta not available, run some adapations first and it will be inferred atomatically")
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

    def adapt(self, result, output):
        if self._has_extension(output, "json"):
            with open(output, "w") as f:
                json.dump(result, output, indent=4)
        if self._has_extension(output, "csv"):
            df = self._to_dataframe(result)
            df.write(output)
        if self._has_extension(output, "tsv"):
            df = self._to_dataframe(result)
            df.write(output, delimiter="\t")
        return result
