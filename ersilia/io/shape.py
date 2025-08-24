class InputShapeSingle(object):
    """
    A class representing a single input shape.
    """

    def __init__(self):
        self.name = "Single"


class InputShape(object):
    """
    A class used to determine the input shape.

    Parameters
    ----------
    input_shape : str or None
        The input shape type. Can only be 'single'. If None, defaults to 'Single'.

    Examples
    --------
    .. code-block:: python

        >>> shape = InputShape("list")
        >>> shape.get().name
        'List'
    """

    def __init__(self, input_shape: str = None):
        if input_shape is None:
            self.shape = InputShapeSingle()
        else:
            self.input_shape = input_shape.lower()
            if self.input_shape == "single":
                self.shape = InputShapeSingle()

    def get(self):
        """
        Get the input shape object.

        Returns
        -------
        object
            The input shape object.
        """
        return self.shape
