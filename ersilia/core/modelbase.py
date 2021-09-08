from ..hub.content.slug import Slug

class ModelBase(object):
    """Base class of a Model.

    """
    def __init__(self, model_id_or_slug):
        slugger = Slug()
        if slugger.is_slug(model_id_or_slug):
            self.slug = model_id_or_slug
            self.model_id = slugger.encode(self.slug)
        else:
            self.model_id = model_id_or_slug
            self.slug = slugger.decode(self.model_id)
