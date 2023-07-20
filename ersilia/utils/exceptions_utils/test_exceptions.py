from .exceptions import ErsiliaError


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
        self.message = (
            "The eos/dest/{0}/information.json file does not exist.".format(
                model_id
            )
        )
        self.hints = (
            "Try fetching and serving the model first."
        )
        super().__init__(self.message, self.hints)


class MissingOutputs(ErsiliaError): 
    def __init__(self, model_id):
        self.message = ("There are not as many outputs as there are inputs. They must be the same.\n")
        self.hints = ("Check whether or not the code for the model is skipping a header, or if any of the smiles used are not working correctly.")
        super().__init__(self.message, self.hints)

