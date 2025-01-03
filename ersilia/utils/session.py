import json
import os
import shutil

import psutil

from ..default import (
    CONTAINER_LOGS_TMP_DIR,
    EOS,
    LOGS_DIR,
    MODELS_JSON,
    SESSION_JSON,
    SESSIONS_DIR,
)


def get_current_pid():
    """
    Get the current process ID.

    Returns
    -------
    int
        The current process ID.
    """
    return os.getpid()


def get_parent_pid():
    """
    Get the parent process ID.

    Returns
    -------
    int
        The parent process ID.
    """
    pid = os.getppid()
    return pid


def get_session_uuid():
    """
    Get the session UUID.

    Returns
    -------
    str
        The session UUID.
    """
    # TODO this should not be implemented here ideally, and callers should use the Session interface in ersilia/core/session.py
    with open(os.path.join(get_session_dir(), SESSION_JSON), "r") as f:
        session = json.load(f)
        return session["identifier"]


def create_session_files(session_name):
    """
    Create session directory and necessary files.

    Parameters
    ----------
    session_name : str
        The name of the session.
    """
    session_dir = os.path.join(SESSIONS_DIR, session_name)
    os.makedirs(os.path.join(session_dir, LOGS_DIR), exist_ok=True)
    os.makedirs(os.path.join(session_dir, CONTAINER_LOGS_TMP_DIR), exist_ok=True)


def create_session_dir():
    """
    Create a session directory.
    """
    remove_orphaned_sessions()
    session_name = f"session_{get_parent_pid()}"
    session_dir = os.path.join(SESSIONS_DIR, session_name)
    os.makedirs(session_dir, exist_ok=True)
    create_session_files(session_name)


def get_session_dir():
    """
    Get the session directory.

    Returns
    -------
    str
        The session directory path.
    """
    return os.path.join(SESSIONS_DIR, get_session_id())


def remove_session_dir(session_name):
    """
    Remove a session directory.

    Parameters
    ----------
    session_name : str
        The name of the session.
    """
    session_dir = os.path.join(SESSIONS_DIR, session_name)
    if os.path.exists(session_dir):
        shutil.rmtree(session_dir)


def determine_orphaned_session():
    """
    Determine orphaned sessions.

    Returns
    -------
    list
        A list of orphaned session names.
    """
    # TODO maybe this is slow, look out for performance
    _sessions = []
    sessions = list(
        filter(lambda s: s.startswith("session_"), os.listdir(SESSIONS_DIR))
    )
    if sessions:
        for session in sessions:
            session_pid = int(session.split("_")[1])
            if not psutil.pid_exists(session_pid):
                _sessions.append(session)
    return _sessions


def remove_orphaned_sessions():
    """
    Remove orphaned sessions.
    """
    orphaned_sessions = determine_orphaned_session()
    for session in orphaned_sessions:
        try:
            remove_session_dir(session)
        except PermissionError:  # TODO Is there a better way to handle this?
            pass  # TODO we should at least log this somehow


def get_session_id():
    """
    Get the session ID.

    Returns
    -------
    str
        The session ID.
    """
    return f"session_{get_parent_pid()}"


def register_model_session(model_id, session_dir):
    """
    Register a model with a session.

    Parameters
    ----------
    model_id : str
        The model ID.
    session_dir : str
        The session directory.
    """
    file_path = os.path.join(EOS, MODELS_JSON)

    if not os.path.exists(file_path):
        with open(file_path, "w") as f:
            json.dump({}, f, indent=4)

    with open(file_path, "r") as f:
        models = json.load(f)

    if (
        model_id not in models
    ):  # TODO This would have implications when we try to run the same model across multiple sessions
        models[model_id] = session_dir
        with open(file_path, "w") as f:
            json.dump(models, f, indent=4)


def get_model_session(model_id):
    """
    Get the model session.

    Parameters
    ----------
    model_id : str
        The model ID.

    Returns
    -------
    str
        The session ID.
    """
    file_path = os.path.join(EOS, MODELS_JSON)
    if not os.path.exists(file_path):
        return None
    with open(file_path, "r") as f:
        models = json.load(f)
    return models.get(model_id, None)


def deregister_model_session(model_id):
    """
    Remove a model from a session.

    Parameters
    ----------
    model_id : str
        The model ID.
    """
    file_path = os.path.join(EOS, MODELS_JSON)
    if not os.path.exists(file_path):
        return
    with open(file_path, "r") as f:
        models = json.load(f)
    if model_id in models:
        del models[model_id]
        with open(file_path, "w") as f:
            json.dump(models, f, indent=4)
