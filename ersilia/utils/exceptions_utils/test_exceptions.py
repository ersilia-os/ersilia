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
