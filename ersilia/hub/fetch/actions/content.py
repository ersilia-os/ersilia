import os
import json
from . import BaseAction
from ....db.hubdata.localslugs import SlugDb
from ...content.card import ModelCard
from ....default import CARD_FILE


class CardGetter(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.mc = ModelCard(config_json=config_json)
        self.slugdb = SlugDb(config_json=config_json)

    def get(self):
        self.logger.debug("Getting model card of {0}".format(self.model_id))
        card = self.mc.get(self.model_id, as_json=False)
        slug = card["Slug"]
        model_path = self._model_path(self.model_id)
        card_path = os.path.join(model_path, CARD_FILE)
        with open(card_path, "w") as f:
            json.dump(card, f, indent=4)
        self.logger.debug("Card saved at {0}".format(card_path))
        self.logger.debug("Saving slug {0}".format(slug))
        self.slugdb.insert(self.model_id, slug)
