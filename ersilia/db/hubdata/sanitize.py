import validators

from ... import ErsiliaBase, throw_ersilia_exception
from ...utils.exceptions_utils.hubdata_exceptions import InvalidUrlInAirtableError
from .interfaces import AirtableInterface


# Potentially related issue: https://github.com/ersilia-os/ersilia/issues/1407
class AirtableSanitizer(ErsiliaBase):
    """
    Sanitizes Airtable records by validating hosted URLs.

    Parameters
    ----------
    config_json : dict
        Configuration settings for initializing the sanitizer.

    """

    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.ai = AirtableInterface(config_json=self.config_json)
        self.HOSTED_URL_FIELD = "Host URL"

    @throw_ersilia_exception()
    def check_hosted_urls(self):
        """
        Checks all hosted URLs in Airtable records and raises an error if any URL is invalid.

        Raises
        ------
        InvalidUrlInAirtableError
            If any URL in the Airtable records is invalid.
        """
        for record in self.ai.items_all():
            fields = record["fields"]
            if self.HOSTED_URL_FIELD in fields:
                url = fields[self.HOSTED_URL_FIELD]
                if not validators.url(url):
                    raise InvalidUrlInAirtableError(url)
