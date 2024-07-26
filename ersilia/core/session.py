import os
import csv
import json
import time
import uuid
import shutil
from ..utils.session import get_session_dir

from ..default import SESSIONS_DIR
from .base import ErsiliaBase


class Session(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        session_dir = get_session_dir()
        self.session_file = os.path.join(session_dir, "session.json")

    def current_model_id(self):
        data = self.get()
        if data is None:
            return None
        else:
            return data["model_id"]

    def current_identifier(self):
        data = self.get()
        if data is None:
            return None
        else:
            return data["identifier"]

    def current_service_class(self):
        data = self.get()
        if data is None or data.get("service_class") is None:
            return None
        else:
            return data["service_class"]

    def register_service_class(self, service_class):
        data = self.get()
        data["service_class"] = service_class
        with open(self.session_file, "w") as f:
            json.dump(data, f, indent=4)

    def tracking_status(self):
        data = self.get()
        if data is None:
            return None
        else:
            return data["track_runs"]

    def open(self, model_id, track_runs):
        self.logger.debug("Opening session {0}".format(self.session_file))
        session = {
            "model_id": model_id,
            "timestamp": str(time.time()),
            "identifier": str(uuid.uuid4()),
            "track_runs": track_runs,
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

    def update_total_memory(self, additional_memory):
        """
        Methods to get the the total memory of processess during serving and running.
        """
        data = self.get()
        if data is None:
            data = {}
        current_memory = float(data.get("total memory used by model(MB)", 0))
        new_memory = current_memory + additional_memory
        data["total memory used by model(MB)"] = f"{new_memory:.5f}"
        with open(self.session_file, "w") as f:
            json.dump(data, f, indent=4)

    def update_cpu_time(self, cpu_time):
        """
        Methods to get the total cpu time of processess during serving and running
        """
        data = self.get()
        if data is None:
            data = {}
        current_cpu = float(data.get("CPU time used by model(seconds)", 0))
        new_cpu = current_cpu + cpu_time
        data["CPU time used by model(seconds)"] = f"{new_cpu}"
        with open(self.session_file, "w") as f:
            json.dump(data, f, indent=4)

    def update_peak_memory(self, peak_memory):
        """
        Update peak memory usage in session data if the new peak is higher.
        """
        data = self.get()
        if "peak memory used by model(MiB)" in data:
            stored_peak_memory = float(data["peak memory used by model(MiB)"])
            if peak_memory > stored_peak_memory:
                data["peak memory used by model(MiB)"] = f"{peak_memory:.5f}"
        else:
            data["peak memory used by model(MiB)"] = f"{peak_memory:.5f}"
        with open(self.session_file, "w") as f:
            json.dump(data, f, indent=4)

    def close(self):
        self.logger.debug("Closing session {0}".format(self.session_file))
        if os.path.isfile(self.session_file):
            os.remove(self.session_file)
