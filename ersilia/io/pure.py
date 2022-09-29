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

    def _is_numeric_array(self):
        if self._is_array():
            data = np.array(self.data).ravel().tolist()
            for x in data:
                if not PureDataTyper(x)._is_numeric():
                    return False
            return True
        else:
            return False

    def _is_string_array(self):
        if self._is_array():
            data = np.array(self.data).ravel().tolist()
            for x in data:
                if not PureDataTyper(x)._is_string():
                    return False
            return True
        else:
            return False

    def _is_mixed_array(self):
        if self._is_array():
            has_numeric = False
            has_string = False
            data = np.array(self.data).ravel().tolist()
            for x in data:
                if PureDataTyper(x)._is_numeric():
                    has_numeric = True
                else:
                    has_string = True
            if has_numeric and has_string:
                return True
            else:
                return False
        else:
            return False

    def get_type(self):
        if self._is_string():
            return {"type": "string"}
        if self._is_numeric():
            return {"type": "numeric"}
        if self._is_numeric_array():
            shape = np.array(self.data).shape
            return {"type": "numeric_array", "shape": shape}
        if self._is_string_array():
            shape = np.array(self.data).shape
            return {"type": "string_array", "shape": shape}
        if self._is_mixed_array():
            shape = np.array(self.data).shape
            return {"type": "mixed_array", "shape": shape}
        return None
