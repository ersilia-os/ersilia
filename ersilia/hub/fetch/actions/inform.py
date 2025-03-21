import json
import os

from ....default import INFORMATION_FILE
from ...bundle.repo import ServiceFile
from ...content.information import Information
from . import BaseAction


class ModelInformer(BaseAction):
    """
    Class to inform about the model by writing information to a JSON file. Contains detail
    metadata about the model.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    config_json : dict
        Configuration settings for the model.
    """

    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.information_file = os.path.join(
            self._dest_dir, self.model_id, INFORMATION_FILE
        )

    def _write_information_json(self):
        data = Information(model_id=self.model_id, config_json=self.config_json).get()
        with open(self.information_file, "w") as f:
            json.dump(data, f, indent=4)

    def _add_info_api(self):
        sf = ServiceFile(
            path=os.path.join(self._get_bundle_location(self.model_id), self.model_id)
        )
        if os.path.exists(sf.get_file()):
            sf.add_info_api(information_file=self.information_file)

    def inform(self):
        """
        Write information to a JSON file and add API info for bentoml models.
        """
        self._write_information_json()
        if self._resolve_pack_method_source(self.model_id) == "bentoml":
            self._add_info_api()
