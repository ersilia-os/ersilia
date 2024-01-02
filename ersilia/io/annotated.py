import numpy as np


class AnnotatedDataTyper(object):
    def __init__(self, data, annotated_type, annotated_shape):
        self.data = data
        self.annotated_type = annotated_type
        self.annotated_shape = annotated_shape

    def _is_string(self):
        if self.annotated_type == "String" and self.annotated_shape == "Single":
            return True
        else:
            return False

    def _is_numeric(self):
        if not self.annotated_shape == "Single":
            return False
        if self.annotated_type not in ["Float", "Integer"]:
            return False
        return True

    def _is_numeric_array(self):
        if self.annotated_shape not in ["List", "Flexible List"]:
            return False
        if self.annotated_type not in ["Float", "Integer"]:
            return False
        return True

    def _is_string_array(self):
        if self.annotated_shape not in ["List", "Flexible List"]:
            return False
        if self.annotated_type != "String":
            return False
        return True

    def _is_mixed_array(self):
        if self.annotated_shape not in ["List", "Flexible List"]:
            return False
        if self.annotated_type not in ["Float", "Integer", "String"]:
            return False
        return True

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
