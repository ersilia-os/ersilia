class InputShapeSingle(object):
    def __init__(self):
        self.name = "Single"


class InputShapeList(object):
    def __init__(self):
        self.name = "List"


class InputShapePairOfLists(object):
    def __init__(self):
        self.name = "Pair of Lists"


class InputShape(object):
    def __init__(self, input_shape):
        if input_shape is None:
            self.shape = InputShapeSingle()
        else:
            self.input_shape = input_shape.lower()
            if self.input_shape == "single":
                self.shape = InputShapeSingle()
            if self.input_shape == "list":
                self.shape = InputShapeList()
            if self.input_shape == "pair of lists":
                self.shape = InputShapePairOfLists()

    def get(self):
        return self.shape
