import json
import numpy as np

try:
    from isaura.core.hdf5 import Hdf5ApiExplorer
except:
    Hdf5ApiExplorer = None

from .base import LakeBase
from ..io.dataframe import Dataframe
from ..io.output import DictlistDataframeConverter


class IsauraInterface(LakeBase):
    def __init__(self, model_id, api_name, config_json):
        LakeBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.api_name = api_name
        if Hdf5ApiExplorer is not None:
            self.hdf5 = Hdf5ApiExplorer(model_id=model_id, api_name=api_name)
            self.converter = DictlistDataframeConverter(config_json=config_json)
            self.is_available = True
        else:
            self.is_available = False

    def _dict_to_list(self, d, input):
        result = []
        for k, i in d.items():
            res = input[i]
            res["idx"] = i
            result += [res]
        return result

    def done_todo(self, input):
        keys = [inp["key"] for inp in input]
        result = self.hdf5.check_keys_exist(key_list=keys)
        done = self._dict_to_list(result["available_keys"], input)
        todo = self._dict_to_list(result["unavailable_keys"], input)
        return done, todo

    def read(self, input):
        keys = [inp["key"] for inp in input]
        values = np.array([res for res in self.hdf5.read_by_key(keys)])
        features = self.hdf5.get_features()
        df = Dataframe(keys=keys, values=np.array(values), features=features)
        results_ = self.converter.dataframe2dictlist(
            df, model_id=self.model_id, api_name=self.api_name
        )
        results = []
        for inp, res in zip(input, results_):
            results += [{"input": inp, "output": res["output"]}]
        return results

    def write(self, results):
        results = json.dumps(results)
        df = self.converter.dictlist2dataframe(
            results, model_id=self.model_id, api_name=self.api_name
        )
        self.hdf5.write_api(df.keys, df.values)
        self.hdf5.write_features(df.features)
