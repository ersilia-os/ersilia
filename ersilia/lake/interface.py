from isaura.core.hdf5 import Hdf5ApiExplorer
from .base import LakeBase


class IsauraInterface(LakeBase):
    def __init__(self, model_id, api_name, config_json):
        LakeBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        self.api_name = api_name
        self.hdf5 = Hdf5ApiExplorer(model_id=model_id, api_name=api_name)

    def _dict_to_lists(self, d):
        keys = []
        vals = []
        for k, v in d.items():
            keys += [k]
            vals += [v]
        return keys, vals

    def done_todo(self, input):
        keys = [inp["key"] for inp in input]
        result = self.hdf5.check_keys_exist(api_name=self.hdf5.api_name, key_list=keys)
        keys, idxs = self._dict_to_lists(result["available_keys"])
        done = {"key": keys, "idx": idxs}
        keys, vals = self._dict_to_lists(result["unavailable_keys"])
        todo = {"key": keys, "idx": idxs}
        return done, todo

    def write(self, results):
        pass

    def read(self, input):
        pass
