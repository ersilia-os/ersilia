import contextlib
import json
import os
import shutil
import stat

try:
    import fcntl
except ImportError:  # not available on Windows
    fcntl = None

import psutil

from ..default import (
    EOS,
    LOGS_DIR,
    MODELS_JSON,
    SESSION_JSON,
    SESSIONS_DIR,
)

os.umask(0)


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
    session_dir = os.path.join(SESSIONS_DIR, session_name)
    logs_dir = os.path.join(session_dir, LOGS_DIR)
    os.makedirs(logs_dir, mode=0o777, exist_ok=True)


def create_session_dir():
    """
    Create a session directory.
    """
    remove_orphaned_sessions()
    session_name = f"session_{get_parent_pid()}"
    session_dir = os.path.join(SESSIONS_DIR, session_name)
    os.makedirs(session_dir, mode=0o777, exist_ok=True)
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


def set_write_permissions(directory):
    current_uid = os.getuid()
    for root, dirs, files in os.walk(directory):
        if os.path.basename(root) == "_logs":
            dirs[:] = []
            continue

        dirs[:] = [d for d in dirs if d != "_logs"]

        for d in dirs:
            dir_path = os.path.join(root, d)
            try:
                if os.stat(dir_path).st_uid == current_uid:
                    os.chmod(dir_path, stat.S_IRWXU)
            except (FileNotFoundError, PermissionError):
                continue

        for f in files:
            file_path = os.path.join(root, f)
            try:
                if os.stat(file_path).st_uid == current_uid:
                    os.chmod(file_path, stat.S_IRWXU)
            except (FileNotFoundError, PermissionError):
                continue


def remove_session_dir(session_name):
    session_dir = os.path.join(SESSIONS_DIR, session_name)
    if os.path.exists(session_dir):
        set_write_permissions(session_dir)
        for item in os.listdir(session_dir):
            item_path = os.path.join(session_dir, item)
            if os.path.basename(item_path) == "_logs":
                continue
            try:
                if os.path.isfile(item_path) or os.path.islink(item_path):
                    os.unlink(item_path)
                elif os.path.isdir(item_path):
                    shutil.rmtree(item_path)
            except FileNotFoundError:
                # Another terminal may be cleaning the same session.
                continue
            except Exception as e:
                raise ValueError(f"Error deleting {item_path}: {e}")
        try:
            shutil.rmtree(session_dir)
        except FileNotFoundError:
            pass
        except Exception as e:
            raise ValueError(f"Error deleting {session_dir}: {e}")


def prune_empty_session_dirs():
    """
    Prune session directories that seem to be empty, meaning they don't contain any serving data or logs.
    """
    for session_name in os.listdir(SESSIONS_DIR):
        session_dir = os.path.join(SESSIONS_DIR, session_name)
        try:
            files = os.listdir(session_dir)
        except (FileNotFoundError, NotADirectoryError):
            continue
        do_prune = True
        for fn in files:
            if fn.endswith(".pid"):
                do_prune = False
                break
            if fn.endswith(".log"):
                do_prune = False
                break
        if do_prune:
            remove_session_dir(session_name)


def session_pid_from_name(session_name):
    """
    Get the parent process ID encoded in a session name.

    Parameters
    ----------
    session_name : str
        A session name such as ``session_12345``, or a path ending in one.

    Returns
    -------
    int or None
        The process ID, or None if the name is not a session name.
    """
    name = os.path.basename(os.path.normpath(session_name))
    prefix, _, pid = name.partition("_")
    if prefix != "session" or not pid.isdigit():
        return None
    return int(pid)


def is_session_alive(session_name):
    """
    Check whether the process that owns a session is still running.

    Parameters
    ----------
    session_name : str
        A session name such as ``session_12345``, or a path ending in one.

    Returns
    -------
    bool
        True if the session's process exists.
    """
    pid = session_pid_from_name(session_name)
    return pid is not None and psutil.pid_exists(pid)


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
            session_pid = session_pid_from_name(session)
            if session_pid is None:
                continue
            if not psutil.pid_exists(session_pid):
                _sessions.append(session)
    return _sessions


def stop_containers_by_name(names):
    """
    Stop and remove the named Docker containers, if any.

    Silent no-op if docker is unavailable or a name is unknown.
    """
    if not names:
        return
    try:
        import docker

        client = docker.from_env()
    except Exception:
        return
    by_name = {c.name: c for c in client.containers.list(all=True)}
    for name in names:
        c = by_name.get(name)
        if c is None:
            continue
        try:
            c.stop()
        except Exception:
            pass
        try:
            c.remove()
        except Exception:
            pass


def kill_process_tree(pid, timeout=5):
    """
    Terminate a process and all of its descendants.

    The model server is launched as ``ersilia_model_serve``, which in turn
    spawns ``run_uvicorn.py`` as a child. Killing only the recorded PID
    orphans the uvicorn server, so the whole tree is terminated instead:
    SIGTERM first, then SIGKILL for anything still alive after ``timeout``.

    Parameters
    ----------
    pid : int
        The root process ID.
    timeout : float, optional
        Seconds to wait for graceful termination before force-killing.
    """
    if pid is None or pid == -1:
        return
    try:
        parent = psutil.Process(pid)
    except psutil.NoSuchProcess:
        return
    try:
        procs = parent.children(recursive=True)
    except psutil.NoSuchProcess:
        procs = []
    procs.append(parent)
    for p in procs:
        try:
            p.terminate()
        except psutil.NoSuchProcess:
            pass
    _, alive = psutil.wait_procs(procs, timeout=timeout)
    for p in alive:
        try:
            p.kill()
        except psutil.NoSuchProcess:
            pass
    if alive:
        psutil.wait_procs(alive, timeout=timeout)


def purge_session_processes(session_name):
    """
    Stop any leftover processes and Docker containers tracked by a session's
    .pid files. Used to clean up after orphaned (terminal-killed) sessions.
    """
    session_dir = os.path.join(SESSIONS_DIR, session_name)
    if not os.path.isdir(session_dir):
        return
    pids = []
    container_names = []
    try:
        files = os.listdir(session_dir)
    except FileNotFoundError:
        return
    for fn in files:
        if not fn.endswith(".pid"):
            continue
        path = os.path.join(session_dir, fn)
        try:
            written_at = os.path.getmtime(path)
            with open(path, "r") as f:
                for line in f:
                    parts = line.strip().split()
                    if not parts:
                        continue
                    try:
                        pids.append((int(parts[0]), written_at))
                    except ValueError:
                        pass
                    if len(parts) >= 3 and parts[2] != "-":
                        container_names.append(parts[2])
        except Exception:
            continue
    for pid, written_at in pids:
        if _pid_was_recycled(pid, written_at):
            continue
        try:
            kill_process_tree(pid)
        except Exception:
            pass
    stop_containers_by_name(container_names)


def _pid_was_recycled(pid, written_at):
    # The recorded server started before its pid file was written. A process
    # with the same pid that started later is an unrelated process (possibly
    # another terminal's model server) and must not be killed.
    if pid is None or pid < 0:
        return False
    try:
        return psutil.Process(pid).create_time() > written_at + 1
    except (psutil.NoSuchProcess, psutil.AccessDenied, ValueError):
        return False


def remove_orphaned_sessions():
    """
    Remove orphaned sessions.
    """
    try:
        orphaned_sessions = determine_orphaned_session()
    except FileNotFoundError:
        return
    for session in orphaned_sessions:
        try:
            purge_session_processes(session)
        except Exception:
            pass
        try:
            remove_session_dir(session)
        except Exception:
            # Cleanup must never stop the CLI, e.g. when another terminal
            # removes the same orphaned session at the same time.
            pass
        try:
            deregister_session(os.path.join(SESSIONS_DIR, session))
        except Exception:
            pass


def get_session_id():
    """
    Get the session ID.

    Returns
    -------
    str
        The session ID.
    """
    return f"session_{get_parent_pid()}"


@contextlib.contextmanager
def _models_json():
    # Yields the {model_id: [session_dir, ...]} mapping and writes it back.
    # A lock serializes the read-modify-write across terminals.
    file_path = os.path.join(EOS, MODELS_JSON)
    lock_path = file_path + ".lock"
    with open(lock_path, "a") as lock:
        if fcntl is not None:
            fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            models = _read_models_json(file_path)
            before = json.dumps(models, sort_keys=True)
            yield models
            models = {k: v for k, v in models.items() if v}
            if json.dumps(models, sort_keys=True) != before:
                tmp_path = file_path + ".tmp"
                with open(tmp_path, "w") as f:
                    json.dump(models, f, indent=4)
                os.replace(tmp_path, file_path)
        finally:
            if fcntl is not None:
                fcntl.flock(lock, fcntl.LOCK_UN)


def _read_models_json(file_path):
    if not os.path.exists(file_path):
        return {}
    try:
        with open(file_path, "r") as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError):
        return {}
    # Older versions stored a single session directory per model.
    return {k: ([v] if isinstance(v, str) else list(v)) for k, v in data.items() if v}


def register_model_session(model_id, session_dir):
    """
    Register that a session is serving a model.

    Several sessions (terminals) can serve the same model at the same time,
    so each model maps to a list of session directories.

    Parameters
    ----------
    model_id : str
        The model ID.
    session_dir : str
        The session directory.
    """
    with _models_json() as models:
        sessions = models.setdefault(model_id, [])
        if session_dir not in sessions:
            sessions.append(session_dir)


def get_model_sessions(model_id):
    """
    Get the sessions registered as serving a model.

    Parameters
    ----------
    model_id : str
        The model ID.

    Returns
    -------
    list
        The session directories, possibly including stale ones.
    """
    file_path = os.path.join(EOS, MODELS_JSON)
    return _read_models_json(file_path).get(model_id, [])


def get_live_model_sessions(model_id):
    """
    Get the sessions that are currently serving a model.

    A session counts as live when its process is running and it still
    holds the model's pid file.

    Parameters
    ----------
    model_id : str
        The model ID.

    Returns
    -------
    list
        The live session directories.
    """
    return [
        session_dir
        for session_dir in get_model_sessions(model_id)
        if is_session_alive(session_dir)
        and os.path.exists(os.path.join(session_dir, f"{model_id}.pid"))
    ]


def deregister_model_session(model_id, session_dir=None):
    """
    Deregister a session from a model, leaving other sessions untouched.

    Parameters
    ----------
    model_id : str
        The model ID.
    session_dir : str, optional
        The session directory to deregister. Defaults to the current session.
    """
    if session_dir is None:
        session_dir = get_session_dir()
    with _models_json() as models:
        sessions = models.get(model_id, [])
        if session_dir in sessions:
            sessions.remove(session_dir)


def deregister_session(session_dir):
    """
    Deregister a session from every model it was registered for.

    Parameters
    ----------
    session_dir : str
        The session directory.
    """
    if not os.path.exists(os.path.join(EOS, MODELS_JSON)):
        return
    with _models_json() as models:
        for sessions in models.values():
            if session_dir in sessions:
                sessions.remove(session_dir)
