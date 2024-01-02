import validators

from .interfaces import AirtableInterface
from ...utils.exceptions_utils.hubdata_exceptions import InvalidUrlInAirtableError

from ... import ErsiliaBase
from ... import throw_ersilia_exception


class AirtableSanitizer(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.ai = AirtableInterface(config_json=self.config_json)
        self.HOSTED_URL_FIELD = "Host URL"

    @throw_ersilia_exception
    def check_hosted_urls(self):
        for record in self.ai.items_all():
            fields = record["fields"]
            if self.HOSTED_URL_FIELD in fields:
                url = fields[self.HOSTED_URL_FIELD]
                if not validators.url(url):
                    raise InvalidUrlInAirtableError(url)
