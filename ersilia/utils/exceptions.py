class ErsiliaError(Exception):
    """Base class for managing errors in Ersilia"""

    def __init__(self, message="Ersilia has experienced an error", hints=""):
        text = "Something went wrong with Ersilia...\n\n"
        text += "{}\n\n".format(self.__class__.__name__)
        text += "Error message:\n"
        text += message
        text += "\n\n"
        if hints:
            text += "Hints:\n"
            text += hints
            text += "\n\n"
        text += "If this error message is not helpful, open an issue at:\n"
        text += " - https://github.com/ersilia-os/ersilia\n"
        text += "Or feel free to reach out to us at:\n"
        text += " - hello[at]ersilia.io\n\n"
        text += (
            "If you haven't, try to run your command in verbose mode (-v in the CLI)"
        )
        super().__init__(text)


class MissingDependencyError(ErsiliaError):
    def __init__(self, dependency):
        self.dependency = dependency
        self.message = "Missing dependency {0}".format(self.dependency)
        self.hints = ""
        super().__init__(self.message, self.hints)


class InvalidModelIdentifierError(ErsiliaError):
    def __init__(self, model):
        self.model = model
        self.message = "Could not identifiy model identifier or slug {0}:".format(
            self.model
        )
        self.hints = "Please check that {0} exists in the Ersilia Model Hub:\n - https://ersilia.io/model-hub (for approved models)\n - https://airtable.com/shrUcrUnd7jB9ChZV (for approved and in preparation models)".format(
            self.model
        )
        super().__init__(self.message, self.hints)


class ModelNotAvailableLocallyError(ErsiliaError):
    def __init__(self, model):
        self.model = model
        self.message = "Model {0} is not available locally, so it cannot be served".format(
            self.model
        )
        self.hints = "Fetch the model using the CLI. Simply run:\n"
        self.hints += "$ ersilia fetch {0}".format(self.model)
        super().__init__(self.message, self.hint)


class EmptyOutputError(ErsiliaError):
    def __init__(self, model_id, api_name):
        self.model_id = model_id
        self.api_name = api_name
        self.message = "Model API {0}:{1} did not produce an output".format(
            self.model_id, self.api_name
        )
        self.hints = "- Visit the fetch troubleshooting site"
        super().__init__(self.message, self.hints)
