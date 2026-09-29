import subprocess
import sys
import time
from unittest.mock import MagicMock, patch

import psutil
import pytest

from ersilia.hub.fetch.fetch import ModelFetcher
from ersilia.serve.autoservice import AutoService
from ersilia.serve.services import _FastApiService
from ersilia.utils.session import kill_process_tree, purge_session_processes

pytestmark = pytest.mark.skipif(
    sys.platform == "win32", reason="POSIX process semantics"
)

# Mirrors ersilia_pack.server:main, which runs run_uvicorn.py via subprocess.run
WRAPPER = "import subprocess, sys; subprocess.run([sys.executable, '-c', {child!r}])"
CHILD_SLEEP = "import time; time.sleep(600)"
CHILD_IGNORES_SIGTERM = (
    "import signal, time; signal.signal(signal.SIGTERM, signal.SIG_IGN); "
    "time.sleep(600)"
)


def _is_running(pid):
    try:
        return psutil.Process(pid).status() != psutil.STATUS_ZOMBIE
    except psutil.NoSuchProcess:
        return False


@pytest.fixture
def server_tree():
    """Spawn a wrapper process with a long-lived child, like ersilia_model_serve."""
    procs = []

    def _spawn(child=CHILD_SLEEP):
        wrapper = subprocess.Popen([sys.executable, "-c", WRAPPER.format(child=child)])
        for _ in range(100):
            children = psutil.Process(wrapper.pid).children()
            if children:
                procs.append((wrapper, children[0].pid))
                return wrapper, children[0].pid
            time.sleep(0.05)
        wrapper.kill()
        wrapper.wait()
        raise RuntimeError("Child process did not start")

    yield _spawn

    # Safety net so a failing test doesn't itself leak processes
    for wrapper, child_pid in procs:
        for pid in (child_pid, wrapper.pid):
            try:
                psutil.Process(pid).kill()
            except psutil.NoSuchProcess:
                pass
        wrapper.wait()


def _assert_tree_gone(wrapper, child_pid):
    wrapper.wait(timeout=10)
    assert not _is_running(wrapper.pid)
    assert not _is_running(child_pid)


def test_kill_process_tree_kills_child(server_tree):
    """Killing the wrapper PID must also kill the server it spawned."""
    wrapper, child_pid = server_tree()
    kill_process_tree(wrapper.pid)
    _assert_tree_gone(wrapper, child_pid)


def test_kill_process_tree_force_kills_after_timeout(server_tree):
    """A child that ignores SIGTERM is force-killed after the timeout."""
    wrapper, child_pid = server_tree(child=CHILD_IGNORES_SIGTERM)
    time.sleep(0.5)  # let the child install its SIGTERM handler
    kill_process_tree(wrapper.pid, timeout=1)
    _assert_tree_gone(wrapper, child_pid)


@pytest.mark.parametrize("pid", [None, -1, 2**22 + 12345])
def test_kill_process_tree_ignores_missing_pids(pid):
    """Unassigned or placeholder PIDs are a no-op, not an error."""
    kill_process_tree(pid)


def test_fastapi_service_close_kills_child(server_tree):
    """_FastApiService.close (used by ersilia close) leaves no orphaned server."""
    wrapper, child_pid = server_tree()
    service = _FastApiService.__new__(_FastApiService)
    service.pid = wrapper.pid
    service.logger = MagicMock()
    service.close()
    _assert_tree_gone(wrapper, child_pid)


def test_autoservice_kill_pids_kills_child(server_tree):
    """AutoService._kill_pids (close / clean_before_serving) kills the whole tree."""
    wrapper, child_pid = server_tree()
    autoservice = AutoService.__new__(AutoService)
    autoservice.logger = MagicMock()
    autoservice._kill_pids([-1, wrapper.pid])
    _assert_tree_gone(wrapper, child_pid)


def test_purge_session_processes_kills_child(server_tree, tmp_path):
    """Orphaned-session cleanup kills the whole tree for PIDs in .pid files."""
    wrapper, child_pid = server_tree()
    session_dir = tmp_path / "session_1"
    session_dir.mkdir()
    (session_dir / "eos0xxx.pid").write_text(
        "{0} http://0.0.0.0:8000 -\n".format(wrapper.pid)
    )
    with patch("ersilia.utils.session.SESSIONS_DIR", str(tmp_path)):
        purge_session_processes("session_1")
    _assert_tree_gone(wrapper, child_pid)


def test_standard_example_closes_model_on_failure():
    """If the fetch-time example run fails, the served model is still closed."""
    fetcher = ModelFetcher.__new__(ModelFetcher)
    fetcher.config_json = None
    example = MagicMock()
    example.run.side_effect = RuntimeError("model crashed")
    with (
        patch("ersilia.hub.fetch.fetch.ModelStandardExample", return_value=example),
        patch("ersilia.hub.fetch.fetch.spinner", side_effect=lambda _, f, **kw: f()),
    ):
        with pytest.raises(RuntimeError):
            fetcher._standard_csv_example("eos0xxx")
    example.close_model.assert_called_once()


def _children_of(pid, timeout=5):
    import psutil

    deadline = time.time() + timeout
    while time.time() < deadline:
        kids = psutil.Process(pid).children()
        if len(kids) >= 2:
            return kids
        time.sleep(0.1)
    raise AssertionError("the helper processes did not start")


def test_a_server_in_its_own_session_is_stopped_with_its_starters():
    # Like Apptainer: the server (the listening process) is not a descendant
    # of the processes that started it, but they share its session's group.
    import psutil

    starter = subprocess.Popen(
        ["sh", "-c", "sleep 60 & sleep 60 & wait"], start_new_session=True
    )
    try:
        server, sibling = _children_of(starter.pid)
        kill_process_tree(server.pid)
        starter.wait(timeout=10)
        assert not psutil.pid_exists(server.pid) or server.status() == "zombie"
        assert not sibling.is_running() or sibling.status() == "zombie"
    finally:
        if starter.poll() is None:
            starter.kill()


def test_a_server_in_our_own_group_leaves_its_siblings_alone():
    # Servers not started in a session of their own: only their own tree.
    first = subprocess.Popen(["sleep", "60"])
    second = subprocess.Popen(["sleep", "60"])
    try:
        kill_process_tree(first.pid)
        first.wait(timeout=10)
        assert second.poll() is None
    finally:
        for p in (first, second):
            if p.poll() is None:
                p.kill()
