import json
import os
import time
import uuid

import psutil

from ..default import SESSION_JSON
from ..utils.session import (
    _pid_was_recycled,
    container_is_running,
    get_session_dir,
    prune_empty_session_dirs,
    read_pid_file,
)
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
        prune_empty_session_dirs()
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

    def current_store_status(self):
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
            return (
                data["read_store"],
                data["write_store"],
                data["store_access"],
                data["nearest_neighbors"],
                data["local_cache"],
            )

    def current_local_cache_status(self):
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
            return data["local_cache"]

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

    def register_store_status(
        self, read_store, write_store, store_access, nearest_neighbors, enable_cache
    ):
        """
        Register the output source in the session.

        This method updates the session data with the provided output source.

        Parameters
        ----------
        output_source : str
            The output source to register.
        """
        data = {} if self.get() is None else self.get()
        data["read_store"] = read_store
        data["write_store"] = write_store
        data["store_access"] = store_access
        data["nearest_neighbors"] = nearest_neighbors
        data["local_cache"] = enable_cache
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
        data = self.get()
        session = {
            "model_id": model_id,
            "timestamp": str(time.time()),
            "identifier": str(uuid.uuid4()),
            "track_runs": track_runs,
        }
        session = session | data if data is not None else session
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

    def served_model(self):
        """
        Tell which model this session serves, and whether it is really running.

        ``session.json`` records which model is served, and ``<model>.pid``
        where it runs. The two can disagree, e.g. when serving was interrupted
        or the model's container was removed outside Ersilia.

        Returns
        -------
        tuple of (str or None, str or None)
            The model ID and its status: "running", or "stale" when the model
            is recorded but its .pid file, container or process is gone.
            (None, None) when no model is recorded.
        """
        data = self.get() or {}
        model_id = data.get("model_id")
        if not model_id:
            return None, None
        pid_file = os.path.join(self._session_dir, "{0}.pid".format(model_id))
        if not os.path.isfile(pid_file):
            return model_id, "stale"
        try:
            pids, containers = read_pid_file(pid_file)
            written_at = os.path.getmtime(pid_file)
        except OSError:
            return model_id, "stale"
        if containers:
            running = [container_is_running(name) for name in containers]
            # If Docker cannot be reached, the record is trusted.
            if any(r is None for r in running) or any(running):
                return model_id, "running"
            return model_id, "stale"
        servers = [p for p in pids if p is not None and p >= 0]
        if servers:
            alive = [
                psutil.pid_exists(p) and not _pid_was_recycled(p, written_at)
                for p in servers
            ]
            return model_id, "running" if any(alive) else "stale"
        # Nothing to check (e.g. a hosted model): trust the record.
        return model_id, "running"

    def clear_stale(self, model_id):
        """
        Forget a model that is recorded as served but is no longer running.

        Parameters
        ----------
        model_id : str
            The model recorded in this session.
        """
        from ..utils.session import deregister_model_session

        self.logger.info("Clearing stale session record of model {0}".format(model_id))
        pid_file = os.path.join(self._session_dir, "{0}.pid".format(model_id))
        if os.path.isfile(pid_file):
            os.remove(pid_file)
        self.close()
        deregister_model_session(model_id, self._session_dir)

    def close(self):
        """
        Close the current session.

        This method removes the session file, effectively closing the session.
        """
        self.logger.debug("Closing session {0}".format(self.session_file))
        if os.path.isfile(self.session_file):
            os.remove(self.session_file)
