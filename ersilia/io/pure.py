import numpy as np

from .. import ErsiliaBase
from ..utils.paths import get_metadata_from_base_dir


class PureDataTyper(ErsiliaBase):
    """
    A class used to determine the type of data provided.

    Parameters
    ----------
    data : any
        The data to be typed.
    model_id : str, optional
        The model identifier, by default None.
    config_json : str, optional
        Path to the configuration JSON file, by default None.

    Examples
    --------
    .. code-block:: python

        >>> data_typer = PureDataTyper(data=[1, 2, 3])
        >>> data_typer.get_type()
        {'type': 'numeric_array', 'shape': (3,)}
    """

    def __init__(self, data: any, model_id: str = None, config_json: str = None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.data = data
        self.model_id = model_id

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
            data = [x for x in data if x is not None]
            if len(data) < 1:
                return False
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
            data = [x for x in data if x is not None]
            if len(data) < 1:
                return False
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

    def get_type_from_metadata(self) -> dict:
        """
        Get the type of data from the model metadata.

        Returns
        -------
        dict
            A dictionary containing the type and shape of the data if available, otherwise None.
        """
        if self.model_id is None:
            return
        dest = self._model_path(self.model_id)
        try:
            meta = get_metadata_from_base_dir(dest)
        except FileNotFoundError:
            return
        output_type = meta["Output Type"]
        output_shape = meta["Output Shape"]
        if len(output_type) > 1:
            return
        if output_shape == "Flexible List":
            return
        output_type = output_type[0]
        if output_shape == "Single":
            if output_type == "Integer":
                return {"type": "numeric"}
            if output_type == "Float":
                return {"type": "numeric"}
            if output_type == "String":
                return {"type": "string"}
            return
        if output_shape == "List":
            if output_type == "Integer":
                return {"type": "numeric_array", "shape": np.array(self.data).shape}
            if output_type == "Float":
                return {"type": "numeric_array", "shape": np.array(self.data).shape}
            if output_type == "String":
                return {"type": "string_array", "shape": np.array(self.data).shape}
            return
        return

    def get_type(self) -> dict:
        """
        Determine the type of the data.

        Returns
        -------
        dict
            A dictionary containing the type and shape of the data.
        """
        data_type = self.get_type_from_metadata()
        if data_type is not None:
            return data_type
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
