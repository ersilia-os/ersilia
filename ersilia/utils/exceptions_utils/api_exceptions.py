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
            "- Verify input matches the model's requirements"
        )
        super().__init__(self.message, self.hints)


class ApiSpecifiedOutputError(ErsiliaError):
    def __init__(self):
        self.message = "Specified output is not correct"
        self.hints = "If you don't specify an output, an interable will be created. If you specify a file extension (.json, .tsv, .csv or .h5), a file will be created. Other valid strings include 'dict', 'numpy', 'pandas' and 'json'"
        ErsiliaError.__init__(self, self.message, self.hints)


# Errors raised by the Python API (ersilia.api). Messages match the CLI's.


class InvalidOptionError(ErsiliaError):
    def __init__(self, message, hints=""):
        self.message = message
        self.hints = hints
        ErsiliaError.__init__(self, self.message, self.hints)


class ModelNotServedError(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = "Model {0} is not being served.".format(model_id)
        self.hints = "Call serve() first."
        ErsiliaError.__init__(self, self.message, self.hints)


class SessionBusyError(ErsiliaError):
    def __init__(self, served_model_id, model_id):
        self.served_model_id = served_model_id
        self.message = "Model {0} is already being served in this session.".format(
            served_model_id
        )
        self.hints = (
            "Close it first, e.g. Model('{0}').close(), then serve {1}.".format(
                served_model_id, model_id
            )
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class ModelFetchError(ErsiliaError):
    def __init__(self, model_id, reason):
        self.model_id = model_id
        self.reason = reason
        self.message = "Model {0} could not be fetched.".format(model_id)
        self.hints = reason or ""
        ErsiliaError.__init__(self, self.message, self.hints)


class ModelServeError(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = "Model {0} could not be started.".format(model_id)
        self.hints = "Try again with verbose=True to see the details."
        ErsiliaError.__init__(self, self.message, self.hints)


class ModelNotDeletableError(ErsiliaError):
    def __init__(self, model_id, reason):
        self.model_id = model_id
        self.message = "Model {0} cannot be deleted.".format(model_id)
        self.hints = reason or ""
        ErsiliaError.__init__(self, self.message, self.hints)


class EmptyRunOutputError(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = "Model {0} produced no output.".format(model_id)
        self.hints = "Check that the input suits this model."
        ErsiliaError.__init__(self, self.message, self.hints)
