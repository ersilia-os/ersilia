"""Errors go to stderr; ERSILIA_SESSION names a session shared by commands."""

import os
import subprocess
import sys

import pytest
from click.testing import CliRunner

import ersilia.utils.echo as echo_utils
import ersilia.utils.session as session_utils
from ersilia.utils.echo import echo, spinner


@pytest.fixture(autouse=True)
def _fresh_echo(monkeypatch):
    monkeypatch.setattr(echo_utils, "_after_error", False)
    monkeypatch.setattr(echo_utils, "_quiet", False)


def _streams(capsys):
    out, err = capsys.readouterr()
    return out, err


# Errors to stderr


def test_an_error_and_its_hints_go_to_stderr(capsys):
    echo("Working on it.")
    echo("Something failed.", fg="red")
    echo("Try this instead.")
    out, err = _streams(capsys)
    assert "Working on it." in out
    assert "Something failed." in err and "Try this instead." in err
    assert "Something failed." not in out and "Try this instead." not in out


def test_output_after_the_error_goes_back_to_stdout(capsys):
    echo("Something failed.", fg="red")
    echo("Model closed.", fg="green")
    echo("Next step.")
    echo("Careful.", fg="yellow")
    out, err = _streams(capsys)
    assert "Something failed." in err
    for line in ("Model closed.", "Next step.", "Careful."):
        assert line in out and line not in err


def test_a_running_step_ends_the_error(capsys):
    echo("Something failed.", fg="red")
    spinner("Working", lambda: None)
    echo("Done.")
    out, err = _streams(capsys)
    assert "Something failed." in err
    assert "Done." in out


def test_a_failed_step_is_printed_to_stderr(capsys):
    def boom():
        raise RuntimeError("no")

    with pytest.raises(RuntimeError):
        spinner("Starting model", boom)
    out, err = _streams(capsys)
    assert "Starting model" in err and "Starting model" not in out


def test_the_cli_writes_errors_to_stderr_only(tmp_path):
    # A real process, so stdout and stderr stay apart.
    env = dict(os.environ, ERSILIA_SESSION="test-streams")
    result = subprocess.run(
        [sys.executable, "-c", "from ersilia.cli import cli; cli()", "info"],
        capture_output=True,
        text=True,
        env=env,
        cwd=tmp_path,
        stdin=subprocess.DEVNULL,
    )
    assert result.returncode == 1
    assert "No model is being served" in result.stderr
    assert "No model is being served" not in result.stdout
    session_utils.remove_session_dir("session_test-streams")


# ERSILIA_SESSION


@pytest.fixture
def eos(tmp_path, monkeypatch):
    sessions = tmp_path / "sessions"
    sessions.mkdir()
    monkeypatch.setattr(session_utils, "EOS", str(tmp_path))
    monkeypatch.setattr(session_utils, "SESSIONS_DIR", str(sessions))
    return tmp_path


def test_the_session_follows_the_terminal_by_default(monkeypatch):
    monkeypatch.delenv("ERSILIA_SESSION", raising=False)
    assert session_utils.get_session_id() == f"session_{os.getppid()}"


def test_ersilia_session_names_the_session(monkeypatch):
    monkeypatch.setenv("ERSILIA_SESSION", "my.project-1")
    assert session_utils.get_session_id() == "session_my.project-1"


@pytest.mark.parametrize("name", ["my project", "12345", "-x", "a/b", "x" * 65])
def test_an_invalid_session_name_is_reported_not_used(monkeypatch, name):
    from ersilia.utils.exceptions_utils.exceptions import SessionNameError

    monkeypatch.setenv("ERSILIA_SESSION", name)
    # While importing, the terminal's session is used; the entry points report it.
    assert session_utils.get_session_id() == f"session_{os.getppid()}"
    assert session_utils.invalid_session_name() == name
    with pytest.raises(SessionNameError, match="not a valid session name"):
        session_utils.check_session_env()


def test_a_named_session_is_alive_and_never_an_orphan(eos):
    named = eos / "sessions" / "session_myproject"
    named.mkdir()
    assert session_utils.is_named_session("session_myproject")
    assert not session_utils.is_named_session("session_12345")
    assert session_utils.is_session_alive(str(named))
    assert "session_myproject" not in session_utils.determine_orphaned_session()


def test_commands_through_a_wrapper_share_a_named_session(tmp_path):
    # Each command below has a different parent process, like 'conda run'.
    env = dict(os.environ, ERSILIA_SESSION="test-shared")
    code = "import ersilia.utils.session as s; print(s.get_session_dir())"
    dirs = {
        subprocess.run(
            ["sh", "-c", f'{sys.executable} -c "{code}"'],
            capture_output=True,
            text=True,
            env=env,
        ).stdout.strip()
        for _ in range(2)
    }
    assert len(dirs) == 1 and dirs.pop().endswith("session_test-shared")


def _serving(eos, name, model_id):
    path = eos / "sessions" / name
    path.mkdir()
    (path / f"{model_id}.pid").write_text("-1 http://0.0.0.0:1 -\n")
    session_utils.register_model_session(model_id, str(path))
    return str(path)


def test_models_served_in_other_sessions_are_listed(eos, monkeypatch):
    here = eos / "sessions" / "session_here"
    here.mkdir()
    monkeypatch.setattr(session_utils, "get_session_dir", lambda: str(here))
    other = _serving(eos, "session_elsewhere", "eos3b5e")
    assert session_utils.models_served_elsewhere() == [("eos3b5e", other)]


def test_no_model_served_points_to_ersilia_session(eos, monkeypatch):
    from ersilia.cli.commands.run import run_cmd

    here = eos / "sessions" / "session_here"
    here.mkdir()
    monkeypatch.setattr(session_utils, "get_session_dir", lambda: str(here))
    _serving(eos, "session_elsewhere", "eos3b5e")
    import ersilia.core.session as core_session

    monkeypatch.setattr(core_session, "get_session_dir", lambda: str(here))
    result = CliRunner().invoke(run_cmd(), ["-i", "in.csv", "-o", "out.csv"])
    assert result.exit_code == 1
    assert "Model eos3b5e is being served in another session" in result.output
    assert "ERSILIA_SESSION=<name>" in result.output


def test_no_hint_when_nothing_is_served_elsewhere(eos, monkeypatch):
    from ersilia.cli.messages import served_elsewhere_hint

    monkeypatch.setattr(session_utils, "get_session_dir", lambda: str(eos))
    served_elsewhere_hint()  # prints nothing and does not fail
