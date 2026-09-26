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


# Stale sessions: session.json names a model that is no longer running.


@pytest.fixture
def session_here(eos, monkeypatch):
    import ersilia.core.session as core_session

    path = _session(eos, os.getpid())
    # A real session always has a log; without one, Session() prunes the dir.
    open(os.path.join(path, "current.log"), "w").close()
    monkeypatch.setattr(core_session, "get_session_dir", lambda: path)
    return path


def _record(path, model_id, pid_line=None):
    with open(os.path.join(path, "session.json"), "w") as f:
        json.dump({"model_id": model_id, "service_class": "pulled_docker"}, f)
    if pid_line is not None:
        with open(os.path.join(path, f"{model_id}.pid"), "w") as f:
            f.write(pid_line + "\n")


def _served(monkeypatch, container_state=True):
    import ersilia.core.session as core_session
    from ersilia.core.session import Session

    monkeypatch.setattr(core_session, "container_is_running", lambda n: container_state)
    return Session(config_json=None).served_model()


def test_served_model_states(session_here, monkeypatch):
    assert _served(monkeypatch) == (None, None)
    _record(session_here, "eos3b5e")
    assert _served(monkeypatch) == ("eos3b5e", "stale")  # no .pid file
    _record(session_here, "eos3b5e", "-1 http://0.0.0.0:1 eos3b5e_1234abcd")
    assert _served(monkeypatch, True) == ("eos3b5e", "running")
    assert _served(monkeypatch, False) == ("eos3b5e", "stale")  # container gone
    assert _served(monkeypatch, None) == ("eos3b5e", "running")  # Docker unknown


def test_served_model_checks_server_processes(session_here, monkeypatch):
    _record(session_here, "eos3b5e", f"{_dead_pid()} http://0.0.0.0:1 -")
    assert _served(monkeypatch) == ("eos3b5e", "stale")
    _record(session_here, "eos3b5e", f"{os.getpid()} http://0.0.0.0:1 -")
    assert _served(monkeypatch) == ("eos3b5e", "running")


def test_clear_stale_forgets_the_model(session_here, monkeypatch):
    from ersilia.core.session import Session

    _record(session_here, "eos3b5e")
    session_utils.register_model_session("eos3b5e", session_here)
    Session(config_json=None).clear_stale("eos3b5e")
    assert not os.path.exists(os.path.join(session_here, "session.json"))
    assert session_utils.get_model_sessions("eos3b5e") == []


@pytest.mark.parametrize(
    "cmd_name, args, expected",
    [
        (
            "run",
            ["-i", "in.csv", "-o", "out.csv"],
            "no longer running in this terminal",
        ),
        ("info", [], "no longer running in this terminal"),
        ("close", [], "closed (it was no longer running)"),
    ],
)
def test_commands_recover_from_a_stale_session(
    session_here, monkeypatch, cmd_name, args, expected
):
    import importlib

    from click.testing import CliRunner

    _record(session_here, "eos3b5e")
    module = importlib.import_module(f"ersilia.cli.commands.{cmd_name}")
    result = CliRunner().invoke(getattr(module, f"{cmd_name}_cmd")(), args)
    assert expected in result.output, result.output
    assert not os.path.exists(os.path.join(session_here, "session.json"))


def test_failed_serve_does_not_leave_a_session_record(monkeypatch):
    from unittest.mock import MagicMock

    import ersilia.core.model as core_model
    from ersilia.core.model import ErsiliaModel

    session = MagicMock()
    monkeypatch.setattr(core_model, "Session", lambda config_json=None: session)
    mdl = ErsiliaModel.__new__(ErsiliaModel)
    mdl.config_json = None
    mdl.model_id = "eos3b5e"
    mdl.logger = logging.getLogger("test")
    mdl.setup = MagicMock()
    mdl.close = MagicMock()
    mdl.autoservice = MagicMock()
    mdl.autoservice.serve.side_effect = KeyboardInterrupt
    with pytest.raises(KeyboardInterrupt):
        ErsiliaModel.serve(mdl)
    session.open.assert_called_once()
    session.close.assert_called_once()
    mdl.autoservice.close.assert_called_once()


def test_serving_the_running_model_again_just_says_so(session_here, monkeypatch):
    from unittest.mock import MagicMock, patch

    from click.testing import CliRunner

    from ersilia.cli.commands.serve import serve_cmd

    _record(session_here, "eos3b5e", "-1 http://0.0.0.0:5000 eos3b5e_1234abcd")
    import ersilia.core.session as core_session

    monkeypatch.setattr(core_session, "container_is_running", lambda name: True)
    monkeypatch.setattr(
        "ersilia.utils.tmp_pid_file",
        lambda model_id: os.path.join(session_here, f"{model_id}.pid"),
    )
    with (
        patch("ersilia.ModelBase", return_value=MagicMock(model_id="eos3b5e")),
        patch("ersilia.core.model.ErsiliaModel") as model,
    ):
        result = CliRunner().invoke(serve_cmd(), ["molecular-weight"])
    assert "Model eos3b5e is already being served in this terminal." in result.output
    assert "It is available at http://127.0.0.1:5000." in result.output
    model.assert_not_called()
