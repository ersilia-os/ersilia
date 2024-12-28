import csv
import os

from ... import ErsiliaBase
from ...default import EOS, FETCHED_MODELS_FILENAME


class FetchedModelsManager(ErsiliaBase):
    """
    Manages fetched models by adding and deleting model entries.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for initializing the fetched models manager.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.file_name = os.path.abspath(os.path.join(EOS, FETCHED_MODELS_FILENAME))

    def add(self, model_id):
        """
        Adds a model entry to the fetched models registry.

        Parameters
        ----------
        model_id : str
            Identifier of the model to add.
        """
        pass

    def delete(self, model_id):
        """
        Deletes a model entry from the fetched models registry(saved in FETCHED_MODELS_FILENAME txt file).

        Parameters
        ----------
        model_id : str
            Identifier of the model to delete.
        """
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
