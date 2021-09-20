import csv
import json
from .. import ErsiliaBase


class DataFrame(object):
    def __init__(self, data, columns):
        self.data = data
        self.columns = columns

    def write(self, file_name, delimiter=None):
        with open(file_name, "w", newline="") as f:
            if delimiter is None:
                writer = csv.writer(f)
            else:
                writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(self.columns)
            for row in self.data:
                writer.writerow(row)


class GenericOutputAdapter(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.logger.debug("Generic output adapter initialized")

    @staticmethod
    def _is_string(output):
        if type(output) is str:
            return True
        else:
            return False

    def _has_extension(self, output, extension):
        if not self._is_string(output):
            return False
        ext = output.split(".")[-1]
        if ext == extension:
            return True
        else:
            return False

    def _to_dataframe(self, result):
        result = json.loads(result)
        R = []
        output_keys = None
        for r in result:
            inp = r["input"]
            out = r["output"]
            if output_keys is None:
                output_keys = [k for k in out.keys()]
            vals = [out[k] for k in output_keys]
            R += [[inp["key"], inp["input"]] + vals]
        df = DataFrame(data=R, columns=["key", "input"] + output_keys)
        return df

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
