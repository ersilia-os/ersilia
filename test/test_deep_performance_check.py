import subprocess
import time
from unittest.mock import patch

import psutil

from ersilia.default import NO_INPUT_PROCESSED_EXIT_CODE
from ersilia.publish.test.services import inspect
from ersilia.publish.test.services.inspect import ModelInspector, _run_with_timeout

RUN_WITH_TIMEOUT = "ersilia.publish.test.services.inspect._run_with_timeout"


def _inspector():
    inspector = ModelInspector.__new__(ModelInspector)
    inspector.model = "eos4qda"
    return inspector


def _completed(returncode):
    return subprocess.CompletedProcess(args="", returncode=returncode)


def _is_gone(pid):
    try:
        return psutil.Process(pid).status() == psutil.STATUS_ZOMBIE
    except psutil.NoSuchProcess:
        return True


def _read_pid(path):
    for _ in range(50):
        if path.exists() and path.read_text().strip():
            return int(path.read_text())
        time.sleep(0.1)
    raise AssertionError("the command did not write its pid")


def test_inputs_that_could_not_be_processed_still_record_the_time():
    with patch(RUN_WITH_TIMEOUT, return_value=NO_INPUT_PROCESSED_EXIT_CODE), patch(
        "subprocess.run", return_value=_completed(0)
    ) as run:
        result, timeout = _inspector()._run_performance_check(1, False)
    assert result.success
    assert "1 predictions executed in" in result.details
    assert not timeout
    # The chain stopped at 'ersilia run', so the model is closed separately.
    assert run.call_args_list[-1].args[0] == "ersilia close"


def test_other_failures_still_fail_the_check():
    with patch(RUN_WITH_TIMEOUT, return_value=1):
        result, timeout = _inspector()._run_performance_check(1, False)
    assert not result.success
    assert not timeout


def test_timeout_records_minus_one_and_closes_the_model():
    with patch(RUN_WITH_TIMEOUT, return_value=None), patch(
        "subprocess.run", return_value=_completed(0)
    ) as run:
        result, timeout = _inspector()._run_performance_check(10, False)
    assert not result.success
    assert result.details == "10 predictions executed in -1.0 seconds. \n"
    assert timeout
    assert run.call_args_list[-1].args[0] == "ersilia close"


def test_after_a_timeout_larger_sizes_are_skipped_without_running():
    with patch(RUN_WITH_TIMEOUT) as run:
        result, timeout = _inspector()._run_performance_check(100, True)
    run.assert_not_called()
    assert result.details == "100 predictions executed in -1.0 seconds. \n"
    assert timeout


def test_run_with_timeout_returns_the_exit_code():
    assert _run_with_timeout("exit 3", timeout=10) == 3


def test_timeout_stops_the_commands_the_shell_started(tmp_path):
    # Like 'ersilia serve && ersilia run': the shell starts a child that
    # would keep running if only the shell were killed.
    pid_file = tmp_path / "child.pid"
    start = time.monotonic()
    code = _run_with_timeout(f"sleep 60 & echo $! > {pid_file}; wait", timeout=1)
    assert code is None
    assert time.monotonic() - start < 10
    assert _is_gone(_read_pid(pid_file))


def test_timeout_kills_a_child_that_ignores_sigterm(tmp_path):
    pid_file = tmp_path / "child.pid"
    start = time.monotonic()
    code = _run_with_timeout(
        f"sh -c 'trap \"\" TERM; echo $$ > {pid_file}; sleep 60' & wait",
        timeout=1,
        grace=1,
    )
    assert code is None
    assert time.monotonic() - start < 10
    assert _is_gone(_read_pid(pid_file))


def test_close_model_does_not_hang_the_check():
    with patch(
        "subprocess.run",
        side_effect=subprocess.TimeoutExpired("ersilia close", 1),
    ):
        inspect._close_model()  # logs a warning instead of raising
