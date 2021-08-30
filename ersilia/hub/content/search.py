"""Search for models"""

from .catalog import CatalogTable


class ModelSearcher(object):
    def __init__(self, catalog):
        self.catalog = catalog

    def search(self, s):
        """Search models based on a string. At the moment, this functionality is very simple."""
        idxs = set()
        s = s.lower()
        for i, r in enumerate(self.catalog.data):
            for r_ in r:
                if s in r_.lower():
                    idxs.update([i])
        R = [r for i, r in enumerate(self.catalog.data) if i in idxs]
        return CatalogTable(data=R, columns=self.catalog.columns)
