from ..shape import InputShapeSingle, InputShapeList, InputShapePairOfLists


class PyInputReader(object):
    def __init__(self, input, IO):
        self.IO = IO
        self.input_shape = IO.input_shape
        if type(self.input_shape) is InputShapeSingle:
            self.expected_number = 1
            self.entity_is_list = False
        if type(self.input_shape) is InputShapeList:
            self.expected_number = 1
            self.entity_is_list = True
        if type(self.input_shape) is InputShapePairOfLists:
            self.expected_number = 2
            self.entity_is_list = True
        if type(input) is tuple:
            input = list(input)
        self._data = input

    def is_single_input(self):
        data = self._data
        if self.entity_is_list:
            assert type(data) is list
            if self.expected_number == 1:
                one_element = data[0]
                if type(one_element) is str:
                    return True
                else:
                    return False
            else:
                one_element = data[0]
                assert type(one_element) is list
                one_inner_element = one_element[0]
                if type(one_inner_element) is list:
                    return False
                else:
                    return True
        else:
            if type(data) is str:
                return True
            else:
                return False

    def read(self):
        if self.is_single_input():
            return [self._data]
        else:
            return self._data
