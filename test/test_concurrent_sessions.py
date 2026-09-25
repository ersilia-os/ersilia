import json
import logging
import os
import subprocess
import sys
import time

import pytest

import ersilia.serve.autoservice as autoservice
import ersilia.utils.session as session_utils
from ersilia.hub.delete.delete import ModelFullDeleter
from ersilia.serve.autoservice import AutoService


@pytest.fixture
def eos(tmp_path, monkeypatch):
    sessions = tmp_path / "sessions"
    sessions.mkdir()
    monkeypatch.setattr(session_utils, "EOS", str(tmp_path))
    monkeypatch.setattr(session_utils, "SESSIONS_DIR", str(sessions))
    return tmp_path


def _dead_pid():
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    proc.wait()
    return proc.pid


def _session(eos, pid, model_id=None):
    path = eos / "sessions" / f"session_{pid}"
    (path / "logs" / "tmp").mkdir(parents=True)
    if model_id:
        (path / f"{model_id}.pid").write_text("-1 http://0.0.0.0:1 -\n")
    return str(path)


def test_clean_temp_dir_only_removes_this_sessions_dirs(eos, tmp_path, monkeypatch):
    mine = _session(eos, 1)
    tmp_root = tmp_path / "tmp"
    own_dir = tmp_root / "ersilia-own"
    other_dir = tmp_root / "ersilia-other"
    own_dir.mkdir(parents=True)
    other_dir.mkdir()
    link = os.path.join(mine, "logs", "tmp", "ersilia-own")
    os.symlink(own_dir, link)
    monkeypatch.setattr(autoservice, "get_session_dir", lambda: mine)

    svc = AutoService.__new__(AutoService)
    svc.logger = logging.getLogger("test")
    svc.clean_temp_dir()

    assert not own_dir.exists()
    assert not os.path.lexists(link)
    assert other_dir.exists()


def test_models_json_keeps_every_session_of_a_model(eos):
    session_utils.register_model_session("eos0aaa", "/s/session_1")
    session_utils.register_model_session("eos0aaa", "/s/session_2")
    session_utils.register_model_session("eos0aaa", "/s/session_2")
    assert session_utils.get_model_sessions("eos0aaa") == [
        "/s/session_1",
        "/s/session_2",
    ]

    session_utils.deregister_model_session("eos0aaa", "/s/session_1")
    assert session_utils.get_model_sessions("eos0aaa") == ["/s/session_2"]

    session_utils.deregister_model_session("eos0aaa", "/s/session_2")
    data = json.loads((eos / "models.json").read_text())
    assert "eos0aaa" not in data


def test_models_json_reads_old_single_session_format(eos):
    (eos / "models.json").write_text(json.dumps({"eos0aaa": "/s/session_1"}))
    assert session_utils.get_model_sessions("eos0aaa") == ["/s/session_1"]
    session_utils.register_model_session("eos0aaa", "/s/session_2")
    assert session_utils.get_model_sessions("eos0aaa") == [
        "/s/session_1",
        "/s/session_2",
    ]


def test_deregister_session_removes_it_from_every_model(eos):
    session_utils.register_model_session("eos0aaa", "/s/session_1")
    session_utils.register_model_session("eos0bbb", "/s/session_1")
    session_utils.register_model_session("eos0bbb", "/s/session_2")
    session_utils.deregister_session("/s/session_1")
    assert session_utils.get_model_sessions("eos0aaa") == []
    assert session_utils.get_model_sessions("eos0bbb") == ["/s/session_2"]


def test_live_model_sessions(eos):
    live = _session(eos, os.getpid(), "eos0aaa")
    dead = _session(eos, _dead_pid(), "eos0aaa")
    closed = _session(eos, os.getppid())
    for s in (live, dead, closed):
        session_utils.register_model_session("eos0aaa", s)
    assert session_utils.get_live_model_sessions("eos0aaa") == [live]


def test_purge_does_not_kill_a_recycled_pid(eos):
    orphan = _session(eos, _dead_pid())
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        pid_file = os.path.join(orphan, "eos0aaa.pid")
        with open(pid_file, "w") as f:
            f.write(f"{proc.pid} http://0.0.0.0:1 -\n")
        # The pid file predates the process, so the pid now belongs to
        # an unrelated process.
        past = time.time() - 600
        os.utime(pid_file, (past, past))

        session_utils.purge_session_processes(os.path.basename(orphan))

        assert proc.poll() is None
    finally:
        proc.kill()
        proc.wait()


def test_purge_kills_the_recorded_server(eos):
    orphan = _session(eos, _dead_pid())
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        # Written after the process started, as serve does.
        time.sleep(1.1)
        with open(os.path.join(orphan, "eos0aaa.pid"), "w") as f:
            f.write(f"{proc.pid} http://0.0.0.0:1 -\n")

        session_utils.purge_session_processes(os.path.basename(orphan))

        assert proc.wait(timeout=10) is not None
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()


def test_orphan_cleanup_tolerates_odd_and_vanished_entries(eos):
    sessions = eos / "sessions"
    (sessions / "session_abc").mkdir()
    (sessions / "session_").mkdir()
    (sessions / "session_123x").write_text("")
    orphan = _session(eos, _dead_pid(), "eos0aaa")
    session_utils.register_model_session("eos0aaa", orphan)

    session_utils.remove_orphaned_sessions()
    session_utils.remove_session_dir("session_999999999")
    session_utils.prune_empty_session_dirs()

    assert not os.path.exists(orphan)
    assert session_utils.get_model_sessions("eos0aaa") == []


def test_delete_refuses_while_another_session_serves_the_model(eos):
    live = _session(eos, os.getpid(), "eos0aaa")
    session_utils.register_model_session("eos0aaa", live)
    md = ModelFullDeleter.__new__(ModelFullDeleter)
    md.logger = logging.getLogger("test")

    can_delete, reason = md.can_be_deleted("eos0aaa")

    assert can_delete is False
    assert "being served" in reason and os.path.basename(live) in reason
    assert os.path.exists(live)


def test_clean_before_serving_stops_containers_of_removed_pid_files(eos, monkeypatch):
    mine = _session(eos, os.getpid())
    with open(os.path.join(mine, "eos0bbb.pid"), "w") as f:
        f.write("-1 http://0.0.0.0:1 eos0bbb_1234abcd\n")
    stopped = []
    monkeypatch.setattr(
        autoservice, "tmp_pid_file", lambda m: os.path.join(mine, f"{m}.pid")
    )
    monkeypatch.setattr(autoservice, "stop_containers_by_name", stopped.extend)

    svc = AutoService.__new__(AutoService)
    svc.logger = logging.getLogger("test")
    svc.model_id = "eos0aaa"
    svc.clean_before_serving()

    assert stopped == ["eos0bbb_1234abcd"]
    assert not os.path.exists(os.path.join(mine, "eos0bbb.pid"))
