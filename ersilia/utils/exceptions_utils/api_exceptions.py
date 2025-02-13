# ruff: noqa: D101, D102
from .exceptions import ErsiliaError


class ApiErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running api command"
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)


class InputFileNotFoundError(ErsiliaError):
    def __init__(self, file_name):
        self.file_name = file_name
        self.message = "Input file {0} does not exist".format(self.file_name)
        self.hints = "Please be make sure that you are passing a valid input file. Accepted formats are .csv, .tsv and .json\n"
        self.hints += "- Check that the file path is correct"
        ErsiliaError.__init__(self, self.message, self.hints)


class UnprocessableInputError(ErsiliaError):
    def __init__(self):
        self.message = "Input data is invalid and cannot be processed"
        self.hints = (
            "No output file will be created.\n"
            "- Check your input data format and content\n"
            "- Ensure chemical structures/identifiers are valid\n"
            "- Verify input matches the model's requirements(valid SMILES)"
        )
        super().__init__(self.message, self.hints)


class ApiSpecifiedOutputError(ErsiliaError):
    def __init__(self):
        self.message = "Specified output is not correct"
        self.hints = "If you don't specify an output, an interable will be created. If you specify a file extension (.json, .tsv, .csv or .h5), a file will be created. Other valid strings include 'dict', 'numpy', 'pandas' and 'json'"
        ErsiliaError.__init__(self, self.message, self.hints)
