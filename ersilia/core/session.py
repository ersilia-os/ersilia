import json
import time
import uuid
import os
import csv
import shutil

from .base import ErsiliaBase
from ..io.output_logger import TabularResultLogger
from ..default import EOS

ERSILIA_RUNS_FOLDER = "ersilia_runs"


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


class RunLogger(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.ersilia_runs_folder = os.path.join(EOS, ERSILIA_RUNS_FOLDER)
        if not os.path.exists(self.ersilia_runs_folder):
            os.mkdir(self.ersilia_runs_folder)
        self.metadata_folder = os.path.join(self.ersilia_runs_folder, "metadata")
        if not os.path.exists(self.metadata_folder):
            os.mkdir(self.metadata_folder)
        self.lake_folder = os.path.join(self.ersilia_runs_folder, "lake")
        if not os.path.exists(self.lake_folder):
            os.mkdir(self.lake_folder)
        self.logs_folder = os.path.join(self.ersilia_runs_folder, "logs")
        if not os.path.exists(self.logs_folder):
            os.mkdir(self.logs_folder)
        self.tabular_result_logger = TabularResultLogger()

    def log_result(self, result):
        output_dir = os.path.join(self.lake_folder, self.model_id)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        file_name = os.path.join(output_dir, "{0}_lake.csv".format(self.model_id))
        tabular_result = self.tabular_result_logger.tabulate(result)
        if tabular_result is None:
            return
        with open(file_name, "w") as f:
            writer = csv.writer(f, delimiter=",")
            for r in tabular_result:
                writer.writerow(r)

    def log_meta(self, meta):
        output_dir = os.path.join(self.metadata_folder, self.model_id)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        file_name = os.path.join(output_dir, "{0}.json".format(self.model_id))
        with open(file_name, "w") as f:
            json.dump(meta, f)

    def log_logs(self):
        output_dir = os.path.join(self.logs_folder, self.model_id)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        file_name = os.path.join(output_dir, "{0}.log".format(self.model_id))
        session_file = os.path.join(EOS, "session.json")
        shutil.copyfile(session_file, file_name)

    def log(self, result, meta):
        self.log_result(result)
        self.log_meta(meta)
        self.log_logs()
