import os
import json

from ...content.information import Information
from ...bundle.repo import ServiceFile
from . import BaseAction

from ....default import INFORMATION_FILE


class ModelInformer(BaseAction):
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
        sf.add_info_api(information_file=self.information_file)

    def inform(self):
        self._write_information_json()
        self._add_info_api()
