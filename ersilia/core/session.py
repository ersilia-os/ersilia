import json
import time
import uuid
import os

from .base import ErsiliaBase
from ..default import EOS


class Session(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.session_file = os.path.join(EOS, "session.json")

    def current_model_id(self):
        data = self.get()
        if data is None:
            return None
        else:
            return data["model_id"]

    def current_service_class(self):
        data = self.get()
        if data is None:
            return None
        else:
            return data["service_class"]

    def register_service_class(self, service_class):
        data = self.get()
        data["service_class"] = service_class
        with open(self.session_file, "w") as f:
            json.dump(data, f, indent=4)

    def open(self, model_id):
        self.logger.debug("Opening session {0}".format(self.session_file))
        session = {
            "model_id": model_id,
            "timestamp": str(time.time()),
            "identifier": str(uuid.uuid4()),
        }
        with open(self.session_file, "w") as f:
            json.dump(session, f, indent=4)

    def get(self):
        if os.path.isfile(self.session_file):
            self.logger.debug("Getting session from {0}".format(self.session_file))
            with open(self.session_file, "r") as f:
                session = json.load(f)
            return session
        else:
            self.logger.debug("No session exists")
            return None

    def close(self):
        self.logger.debug("Closing session {0}".format(self.session_file))
        if os.path.isfile(self.session_file):
            os.remove(self.session_file)
