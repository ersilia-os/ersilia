import numpy as np


class AnnotatedDataTyper(object):
    """
    Class to determine the type of annotated data.

    Parameters
    ----------
    data : any
        The annotated data.
    annotated_type : str
        The type of the annotated data (e.g., "String", "Float", "Integer").
    annotated_shape : str
        The shape of the annotated data (e.g., "Single", "List", "Flexible List").
    """

    def __init__(self, data, annotated_type, annotated_shape):
        self.data = data
        self.annotated_type = annotated_type
        self.annotated_shape = annotated_shape

    def get_type(self):
        """
        Get the type of the annotated data.

        Returns
        -------
        dict or None
            Dictionary containing the type and shape of the annotated data, or None if the type is not recognized.
        """
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
