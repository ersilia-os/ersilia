import socket
from contextlib import closing


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
