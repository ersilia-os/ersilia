from ... import ErsiliaBase
from ...db.hubdata.localslugs import SlugDb
from ...utils.identifiers.model import ModelIdentifier
from .card import ModelCard


class Slug(ErsiliaBase):
    """
    Class to handle slug operations.

    This class provides methods to encode and decode slugs to model IDs
    and vice versa, using local and remote sources.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings in JSON format.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.db = SlugDb(config_json=config_json)
        self.mi = ModelIdentifier()
        self.mc = ModelCard(config_json=config_json)

    def is_slug(self, text: str) -> bool:
        """
        Check if the given text is a slug.

        Parameters
        ----------
        text : str
            The text to check.

        Returns
        -------
        bool
            True if the text is a slug, False otherwise.
        """
        if self.mi.is_valid(text):
            return False
        else:
            return True

    def _local_encode(self, slug):
        r = self.db.models_of_slug(slug)
        if not r:
            self.db.delete_by_slug(slug)
            return None
        else:
            return list(r)[0]

    def _remote_encode(self, slug):
        res = self.mc.get(slug=slug)
        if res is None:
            return None
        else:
            return res["Identifier"].strip()

    def encode(self, slug: str) -> str:
        """
        Given a slug, get the model_id.

        Parameters
        ----------
        slug : str
            The slug to encode.

        Returns
        -------
        str
            The model_id corresponding to the slug.
        """
        model_id = self._local_encode(slug)
        if model_id is not None:
            return model_id
        model_id = self._remote_encode(slug)
        return model_id

    def _local_decode(self, model_id):
        r = self.db.slugs_of_model(model_id)
        if not r:
            self.db.delete_by_model_id(model_id)
            return None
        else:
            return list(r)[0]

    def _remote_decode(self, model_id):
        res = self.mc.get(model_id=model_id)
        if res is None:
            return None
        else:
            if "card" in res:
                res = res["card"]
            return res["Slug"].strip()

    def decode(self, model_id: str) -> str:
        """
        Given a model_id, get the slug.

        Parameters
        ----------
        model_id : str
            The model_id to decode.

        Returns
        -------
        str
            The slug corresponding to the model_id.
        """
        slug = self._local_decode(model_id)
        if slug is not None:
            return slug
        slug = self._remote_decode(model_id)
        return slug
