from .exceptions import ErsiliaError


class InvalidUrlInAirtableError(ErsiliaError):
    def __init__(self, url):
        self.message = "There is an invalid URL in AirTable: {0}".format(url)
        self.hints = "Check the Ersilia Model Hub AirTable database and correct the error. You will need write permissions."
        ErsiliaError.__init__(self, self.message, self.hints)
