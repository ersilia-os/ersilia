import requests
import os
import json
import collections
from .. import ErsiliaBase


class ModelCardSchema(object):

    def __init__(self, model_id, data):
        self.model_id = model_id
        data_ = data[model_id]
        self.title = data_["title"]
        self.date = data_["date"]
        self.author = data_["author"]
        self.categories = data_["categories"]
        self.tags = data_["tags"]
        self.hidden = data_["hidden"]
        self.featured = data_["featured"]
        self.description = data_["description"]
        self.summary = data_["summary"]
        self.github = data_["github"]
        self.app = data_["app"]
        self.publication = data_["publication"]
        odata = collections.OrderedDict()
        odata["model_id"] = self.model_id
        for k in ["title", "date", "author", "categories", "tags",
                  "description", "summary", "app", "github", "publication"]:
            odata[k] = data_[k]
        self._data = odata

    def __str__(self):
        return self.as_json()

    def __repr__(self):
        return self.as_json()

    def as_dict(self):
        return self._data

    def as_json(self):
        return json.dumps(self.as_dict(), indent=4)


class ModelCard(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.url = self.cfg.HUB.WEB
        self.data = json.loads(requests.get(os.path.join(self.url, "db.json")).json())

    def get(self, model_id):
        return ModelCardSchema(model_id, self.data)
