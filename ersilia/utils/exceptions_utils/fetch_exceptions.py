from .exceptions import ErsiliaError

class FetchErsiliaError(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = "Error occured while fetching model: {0}".format(self.model_id)
        self.hints = ""
        super().__init__(self.message, self.hints)

class GetFetchErsiliaError(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = "Error occured while fetching model: {0}".format(self.model_id)
        self.hints = ""
        super().__init__(self.message, self.hints)