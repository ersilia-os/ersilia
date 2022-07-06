import os
import csv

from ... import ErsiliaBase
from ...default import FETCHED_MODELS_FILENAME, EOS


class FetchedModelsManager(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.file_name = os.path.abspath(os.path.join(EOS, FETCHED_MODELS_FILENAME))

    def add(self, model_id):
        pass

    def delete(self, model_id):
        if os.path.exists(self.file_name):
            with open(self.file_name) as infile:
                models = dict(csv.reader(infile))
            infile.close()
            del models[model_id]
            with open(self.file_name, "w") as f:
                for key, values in models.items():
                    f.write(f"{key},{values}\n")
            self.logger.debug("Fetched model entry {0} deleted".format(model_id))
        else:
            self.logger.debug(
                "Model entry {0} was not available in the fetched models registry".format(
                    model_id
                )
            )
