from .exceptions import ErsiliaError


class EmptyInputError(ErsiliaError):
    def __init__(self):
        self.message = "Some items in your input appear to be empty!"
        self.hints = "Please make sure that no empty lines are present in your input. If you are using a file as input and you can't find the empty line, look at the end of the file! It is possible that it contains an extra empty line."
        super().__init__(self.message, self.hints)


class ColumnIndexInputError(ErsiliaError):
    def __init__(self):
        self.message = "The column index could not be accessed!"
        self.hints = "There is some inconsistency in your tabular input. A column could not be accessed. Please check your input table."
        super().__init__(self.message, self.hints)
