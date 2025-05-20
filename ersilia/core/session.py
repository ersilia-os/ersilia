import json
import os
import time
import uuid

from ..default import SESSION_JSON
from ..utils.session import get_session_dir
from .base import ErsiliaBase


class Session(ErsiliaBase):
    """
    Session class for managing model sessions.

    This class provides functionality to manage sessions, including opening, closing,
    and updating session information. Sessions are essential for tracking the state
    and usage of models, ensuring that all necessary information is stored and can be
    retrieved when needed.

    Parameters
    ----------
    config_json : dict
        Configuration in JSON format.
    """

    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self._session_dir = get_session_dir()
        self.session_file = os.path.join(self._session_dir, SESSION_JSON)

    def current_model_id(self):
        """
        Get the current model ID from the session.

        This method retrieves the current model ID from the session data.

        Returns
        -------
        str or None
            The current model ID, or None if no session data is available.
        """
        data = self.get()
        if data is None:
            return None
        else:
            return data["model_id"]

    def current_identifier(self):
        """
        Get the current identifier from the session.

        This method retrieves the current identifier from the session data.

        Returns
        -------
        str or None
            The current identifier, or None if no session data is available.
        """
        data = self.get()
        if data is None:
            return None
        else:
            return data["identifier"]

    def current_service_class(self):
        """
        Get the current service class from the session.

        This method retrieves the current service class from the session data.

        Returns
        -------
        str or None
            The current service class, or None if no session data is available.
        """
        data = self.get()
        if data is None or data.get("service_class") is None:
            return None
        else:
            return data["service_class"]

    def current_output_source(self):
        """
        Get the current output source from the session.

        This method retrieves the current output source from the session data.

        Returns
        -------
        str or None
            The current output source, or None if no session data is available.
        """
        data = self.get()
        if data is None:
            return None
        else:
            return data["output_source"]

    def register_service_class(self, service_class):
        """
        Register the service class in the session.

        This method updates the session data with the provided service class.

        Parameters
        ----------
        service_class : str
            The service class to register.
        """
        data = self.get()
        data["service_class"] = service_class
        with open(self.session_file, "w") as f:
            json.dump(data, f, indent=4)

    def register_output_source(self, output_source):
        """
        Register the output source in the session.

        This method updates the session data with the provided output source.

        Parameters
        ----------
        output_source : str
            The output source to register.
        """
        data = self.get()
        data["output_source"] = output_source
        with open(self.session_file, "w") as f:
            json.dump(data, f, indent=4)

    def register_tracking_use_case(self, use_case):
        """
        Register the tracking use case in the session.

        This method updates the session data with the provided tracking use case.

        Parameters
        ----------
        use_case : str
            The tracking use case to register.
        """
        data = self.get()
        data["tracking_use_case"] = use_case
        with open(self.session_file, "w") as f:
            json.dump(data, f, indent=4)
        self.logger.debug("Registering tracking use case: {0}".format(use_case))

    def tracking_status(self):
        """
        Get the tracking status from the session.

        This method retrieves the tracking status from the session data.

        Returns
        -------
        bool or None
            The tracking status, or None if no session data is available.
        """
        data = self.get()
        if data is None:
            return None
        else:
            return data["track_runs"]

    def get_tracking_use_case(self):
        """
        Get the tracking use case from the session.

        This method retrieves the tracking use case from the session data.

        Returns
        -------
        str or None
            The tracking use case, or None if no session data is available.
        """
        data = self.get()
        if data is None:
            return None
        else:
            return data.get("tracking_use_case", None)

    def open(self, model_id, track_runs):
        """
        Open a new session for the specified model.

        This method creates a new session for the specified model and saves the session data.

        Parameters
        ----------
        model_id : str
            The identifier of the model.
        track_runs : bool
            Whether to track runs.
        """
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
        """
        Get the current session data.

        This method retrieves the current session data from the session file. The session file
        is a JSON file that contains information about the current session, such as the model ID,
        timestamp, identifier, tracking status, service class, and output source.

        Returns
        -------
        dict or None
            The session data, or None if no session file exists.
        """
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
        Update the total memory usage in the session data.

        This method updates the total memory usage in the session data by adding the provided
        additional memory.

        Parameters
        ----------
        additional_memory : float
            The additional memory to add.
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
        Updates the total CPU time usage in the session data by adding the provided
        CPU time.

        Parameters
        ----------
        cpu_time : float
            The CPU time to add.
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
        Update the peak memory usage in the session data.

        This method updates the peak memory usage in the session data if the new peak is higher
        than the stored peak memory.

        Parameters
        ----------
        peak_memory : float
            The peak memory usage to update.
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
        """
        Close the current session.

        This method removes the session file, effectively closing the session.
        """
        self.logger.debug("Closing session {0}".format(self.session_file))
        if os.path.isfile(self.session_file):
            os.remove(self.session_file)
