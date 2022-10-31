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
