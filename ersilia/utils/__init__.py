import os
from ..default import SESSIONS_DIR
from ..utils.session import get_session_id


def tmp_pid_file(model_id):
    """Gets filename to store process ID"""
    session_dir = os.path.join(SESSIONS_DIR, get_session_id())
    return os.path.join(session_dir, "{0}.pid".format(model_id))
