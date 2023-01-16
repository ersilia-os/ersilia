from .exceptions import ErsiliaError


class ApiErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running api command"
        self.hints = ""
        super().__init__(self.message, self.hints)


class InputFileNotFoundError(ErsiliaError):
    def __init__(self, file_name):
        self.file_name = file_name
        self.message = "Input file {0} does not exist".format(self.file_name)
        self.hints = "Please be make sure that you are passing a valid input file. Accepted formats are .csv, .tsv and .json\n"
        self.hints += "- Check that the file path is correct"
        super().__init__(self.message, self.hints)


class ApiSpecifiedOutputError(ErsiliaError):
    def __init__(self):
        self.message = "Specified output is not correct"
        self.hints = "If you don't specify an output, an interable will be created. If you specify a file extension (.json, .tsv, .csv or .h5), a file will be created. Other valid strings include 'dict', 'numpy', 'pandas' and 'json'"
        super().__init__(self.message, self.hints)
