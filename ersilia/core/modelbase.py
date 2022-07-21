import os
import sys
import json

from .. import ErsiliaBase
from ..hub.content.slug import Slug
from ..hub.fetch import STATUS_FILE, DONE_TAG

from ..utils.exceptions import InvalidModelIdentifierError


class ModelBase(ErsiliaBase):
    """Base class of a Model."""

    def __init__(self, model_id_or_slug, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.text = model_id_or_slug
        slugger = Slug()
        try:
            if slugger.is_slug(model_id_or_slug):
                self.slug = model_id_or_slug
                self.model_id = slugger.encode(self.slug)
            else:
                self.model_id = model_id_or_slug
                self.slug = slugger.decode(self.model_id)
            if not self.is_valid():
                raise InvalidModelIdentifierError(model=self.text)
        except InvalidModelIdentifierError as err:
            self.logger.error("Invalid model identifier!")
            print(err)
            sys.exit()
        finally:
            pass

    def is_valid(self):
        if self.model_id is None or self.slug is None:
            return False
        else:
            return True

    def is_available_locally(self):
        fetch_status_file = os.path.join(self._dest_dir, self.model_id, STATUS_FILE)
        if not os.path.exists(fetch_status_file):
            self.logger.debug("No status file exists")
            is_fetched = False
        else:
            with open(fetch_status_file, "r") as f:
                status = json.load(f)
            is_fetched = status[DONE_TAG]
        self.logger.debug("Is fetched: {0}".format(is_fetched))
        return is_fetched
