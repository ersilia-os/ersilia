import os
import shutil
import psutil

from ..default import SESSIONS_DIR, LOGS_DIR, CONTAINER_LOGS_TMP_DIR

def get_current_pid():
    return os.getpid()

def get_parent_pid():
    pid = os.getppid()
    return pid

def create_session_files(session_name):
    session_dir = os.path.join(SESSIONS_DIR, session_name)
    os.makedirs(os.path.join(session_dir, LOGS_DIR), exist_ok=True)
    os.makedirs(os.path.join(session_dir, CONTAINER_LOGS_TMP_DIR), exist_ok=True)

def create_session_dir():
    remove_orphaned_sessions()
    session_name = f"session_{get_parent_pid()}"
    session_dir = os.path.join(SESSIONS_DIR, session_name)
    os.makedirs(session_dir, exist_ok=True)
    create_session_files(session_name)

def get_session_dir():
    return os.path.join(SESSIONS_DIR, get_session_id())

def remove_session_dir(session_name):
    session_dir = os.path.join(SESSIONS_DIR, session_name)
    shutil.rmtree(session_dir)

def determine_orphaned_session():
    # TODO maybe this is slow, look out for performance
    _sessions = []
    sessions = list(filter(lambda s: s.startswith("session_"), os.listdir(SESSIONS_DIR)))
    if sessions:
        for session in sessions:
            session_pid = int(session.split("_")[1])
            if not psutil.pid_exists(session_pid):
                _sessions.append(session)
    return _sessions

def remove_orphaned_sessions():
    orphaned_sessions = determine_orphaned_session()
    for session in orphaned_sessions:
        try:
            remove_session_dir(session)
        except PermissionError:  #TODO Is there a better way to handle this?
            pass  # TODO we should at least log this somehow

def get_session_id():
    return f"session_{get_parent_pid()}"

