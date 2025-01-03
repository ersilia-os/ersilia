import os

from ..utils.session import get_session_dir


def tmp_pid_file(model_id):
    """Gets filename to store process ID"""
    session_dir = get_session_dir()
    return os.path.join(session_dir, "{0}.pid".format(model_id))
