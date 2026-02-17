import socket
import time
from contextlib import closing

import requests


def is_port_in_use(port):
    """
    Check if a port is in use.

    Parameters
    ----------
    port : int
        The port number to check.

    Returns
    -------
    bool
        True if the port is in use, False otherwise.
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(("localhost", port)) == 0


def find_free_port(preferred_port=None):
    """
    Find a free port.

    Parameters
    ----------
    preferred_port : int, optional
        The preferred port number to check first. Default is None.

    Returns
    -------
    int
        A free port number.
    """
    if preferred_port is not None:
        if not is_port_in_use(preferred_port):
            return preferred_port
    with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        s.bind(("", 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]


def _ensure_ready(self, root, attempts=60, sleep_s=0.5):
    if getattr(self, "_ready", False):
        return
    for i in range(1, attempts + 1):
        try:
            r = requests.get(root, timeout=2)
            if 200 <= r.status_code < 500:
                self._ready = True
                self.logger.info(f"Probe OK {root} on attempt {i}")
                return
            self.logger.info(f"Probe {root} attempt {i} got {r.status_code}")
        except Exception as e:
            self.logger.info(f"Probe {root} attempt {i} error: {e}")
        time.sleep(sleep_s)
    raise RuntimeError(f"Server not ready after {attempts} attempts at {root}")
