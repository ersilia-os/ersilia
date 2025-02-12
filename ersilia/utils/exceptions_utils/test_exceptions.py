from .exceptions import ErsiliaError

# ruff: noqa: D101, D102


class WrongCardIdentifierError(ErsiliaError):
    def __init__(self, model_id):
        self.message = (
            "The model identifier in the model card is not correct: {0}".format(
                model_id
            )
        )
        self.hints = (
            "Check the model information, usually available in a metadata.json file."
        )
        super().__init__(self.message, self.hints)


class InformationFileNotExist(ErsiliaError):
    def __init__(self, model_id):
        self.message = "The eos/dest/{0}/information.json file does not exist.".format(
            model_id
        )
        self.hints = "Try fetching and serving the model first, and make sure the model is written correctly."
        super().__init__(self.message, self.hints)


class MissingOutputs(ErsiliaError):
    def __init__(self):
        self.message = "There are not as many outputs as there are inputs. They must be the same.\n"
        self.hints = "Try checking whether or not the code for the model is skipping a header, or if any of the smiles used are not working correctly."
        super().__init__(self.message, self.hints)


class InconsistentOutputs(ErsiliaError):
    def __init__(self, model_id):
        self.message = "Model outputs are inconsistent, meaning the outputs do not fit the 5% similarity threshold."
        self.hints = "Observe the output comparisons for each smiles input above."
        super().__init__(self.message, self.hints)


class EmptyField(ErsiliaError):
    def __init__(self, empty_field):
        self.message = "The {0} field in the model card is empty.".format(empty_field)
        self.hints = (
            "Check the model information, usually available in a metadata.json file."
        )
        super().__init__(self.message, self.hints)


class EmptyKey(ErsiliaError):
    def __init__(self, empty_field):
        self.message = "The {0} key in the model card is empty.".format(empty_field)
        self.hints = "Check the model information, usually available in a metadata.json file or metadata.yml file."
        super().__init__(self.message, self.hints)


class InvalidEntry(ErsiliaError):
    def __init__(self, invalid_field):
        self.message = (
            "The {0} field of this model is not recognized as a valid format.".format(
                invalid_field
            )
        )
        self.hints = (
            "Check the model information, usually available in a metadata.json file."
        )
        super().__init__(self.message, self.hints)


class InconsistentOutputTypes(ErsiliaError):
    def __init__(self, model_id):
        self.message = "Model output types are inconsistent."
        self.hints = "Observe the output comparisons above for each input, and pay attention to the type of the output (string, float, list, etc.) because they do not match."
        super().__init__(self.message, self.hints)
