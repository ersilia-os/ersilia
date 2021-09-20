import numpy as np


class PureDataTyper(object):
    def __init__(self, data):
        self.data = data

    def _is_string(self):
        if type(self.data) is str:
            return True
        else:
            return False

    def _is_numeric(self):
        try:
            float(self.data)
            return True
        except:
            return False

    def _is_array(self):
        try:
            np.array(self.data)
            return True
        except:
            return False

    def get_type(self):
        if self._is_string():
            return {"type": "string"}
        if self._is_numeric():
            return {"type": "numeric"}
        if self._is_array():
            shape = np.array(self.data).shape
            return {"type": "array", "shape": shape}
        return None
