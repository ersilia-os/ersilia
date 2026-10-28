import logging
import os

import pytest

import ersilia.serve.autoservice as autoservice
import ersilia.utils.conda as conda_utils
import ersilia.utils.venv as venv_utils
from ersilia.serve.autoservice import AutoService
from ersilia.utils.conda import SimpleConda
from ersilia.utils.exceptions_utils.fetch_exceptions import ModelPackageInstallError
from ersilia.utils.exceptions_utils.throw_ersilia_exception import (
    is_library_mode,
    set_library_mode,
)
from ersilia.utils.logging import (
    KEEP_TEMP_DIR_MARKER,
    is_kept_temp_dir,
    keep_temp_dir,
)
from ersilia.utils.venv import SimpleVenv

FAILED_INSTALL_LOG = "ERROR: No matching distribution found for x\n"
# pip prints this line even when the install succeeds, so it is not a failure.
RESOLVER_WARNING_LOG = (
    "ERROR: pip's dependency resolver does not currently take into account "
    "all the packages that are installed.\n"
)


@pytest.fixture
def library_mode():
    # run_commandlines is wrapped in throw_ersilia_exception, which prints the
    # error and exits; library mode makes it raise so the test can check it.
    previous = is_library_mode()
    set_library_mode(True)
    yield
    set_library_mode(previous)


def _session_with_tmp_dirs(tmp_path, names):
    session = tmp_path / "session"
    logs_tmp = session / "logs" / "tmp"
    logs_tmp.mkdir(parents=True)
    dirs = {}
    for name in names:
        target = tmp_path / "tmp" / name
        target.mkdir(parents=True)
        os.symlink(target, logs_tmp / name)
        dirs[name] = target
    return str(session), dirs


def _fake_shell(monkeypatch, module, tmp_path, log_text):
    # Use a fixed temp folder and, instead of running the command, write
    # log_text to the log file that the command pipes into with "tee -a".
    folder = tmp_path / "ersilia-run"
    folder.mkdir()
    monkeypatch.setattr(module, "make_temp_dir", lambda prefix: str(folder))

    def fake_run(cmd, quiet=True):
        log = cmd.split("tee -a ")[1].strip()
        with open(log, "w") as f:
            f.write(log_text)

    monkeypatch.setattr(module, "run_command", fake_run)
    return folder


def _conda(monkeypatch):
    conda = SimpleConda.__new__(SimpleConda)
    monkeypatch.setattr(conda, "exists", lambda env: True)
    monkeypatch.setattr(
        conda,
        "create_executable_bash_script",
        lambda env, cmds, fn: open(fn, "w").close(),
    )
    return conda


def test_keep_temp_dir_marks_the_dir(tmp_path):
    assert not is_kept_temp_dir(str(tmp_path))
    keep_temp_dir(str(tmp_path))
    assert (tmp_path / KEEP_TEMP_DIR_MARKER).exists()
    assert is_kept_temp_dir(str(tmp_path))


def test_keep_temp_dir_never_raises(tmp_path):
    missing = tmp_path / "gone"
    keep_temp_dir(str(missing))
    assert not is_kept_temp_dir(str(missing))


def test_clean_temp_dir_keeps_marked_dirs(tmp_path, monkeypatch):
    session, dirs = _session_with_tmp_dirs(tmp_path, ["ersilia-failed", "ersilia-ok"])
    (dirs["ersilia-failed"] / "command_outputs.log").write_text(FAILED_INSTALL_LOG)
    keep_temp_dir(str(dirs["ersilia-failed"]))
    monkeypatch.setattr(autoservice, "get_session_dir", lambda: session)
    svc = AutoService.__new__(AutoService)
    svc.logger = logging.getLogger("test")

    svc.clean_temp_dir()

    assert (dirs["ersilia-failed"] / "command_outputs.log").exists()
    assert not dirs["ersilia-ok"].exists()


def test_install_error_hints_show_the_log_path():
    err = ModelPackageInstallError("bash script.sh", log_file="/tmp/x/out.log")
    assert "/tmp/x/out.log" in err.hints


def test_install_error_without_log_path_keeps_old_hints():
    err = ModelPackageInstallError("bash script.sh")
    assert err.hints.startswith("Try to manually activate")


def test_conda_failure_keeps_log_and_reports_it(tmp_path, monkeypatch, library_mode):
    folder = _fake_shell(monkeypatch, conda_utils, tmp_path, FAILED_INSTALL_LOG)

    with pytest.raises(ModelPackageInstallError) as info:
        _conda(monkeypatch).run_commandlines("eos0xxx", "pip install x")

    log = str(folder / "command_outputs.log")
    assert info.value.log_file == log
    assert log in info.value.hints
    assert is_kept_temp_dir(str(folder))


def test_conda_resolver_warning_does_not_keep_dir(tmp_path, monkeypatch, library_mode):
    folder = _fake_shell(monkeypatch, conda_utils, tmp_path, RESOLVER_WARNING_LOG)

    _conda(monkeypatch).run_commandlines("eos0xxx", "pip install x")

    assert not is_kept_temp_dir(str(folder))


def test_venv_failure_keeps_log_and_reports_it(tmp_path, monkeypatch, library_mode):
    folder = _fake_shell(monkeypatch, venv_utils, tmp_path, FAILED_INSTALL_LOG)
    venv = SimpleVenv.__new__(SimpleVenv)
    venv.root = str(tmp_path)
    monkeypatch.setattr(venv, "exists", lambda env: True)

    with pytest.raises(ModelPackageInstallError) as info:
        venv.run_commandlines("eos0xxx", "pip install x")

    assert info.value.log_file == str(folder / "installs.log")
    assert is_kept_temp_dir(str(folder))
