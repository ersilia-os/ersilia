import os
import re
import stat
import sys

import pytest

from ersilia.publish.test.services import inspect as inspect_module
from ersilia.publish.test.services.inspect import ModelInspector

pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="POSIX shell")

# Mirrors parse_performance in ersilia/publish/test/services/io.py
PERFORMANCE_PATTERN = r"(\d+) predictions executed in (-?\d+\.\d+) seconds. \n"

# Stands in for the ersilia CLI: logs each call, and fails one subcommand on request.
FAKE_ERSILIA = """#!/bin/sh
echo "$1" >> "$FAKE_ERSILIA_LOG"
if [ "$1" = "close" ] && [ "$FAKE_CLOSE_FAIL" = "1" ]; then
    echo "Model eos0abc could not be closed." >&2
    exit 1
fi
if [ "$1" = "$FAKE_FAIL_CMD" ]; then
    echo "$FAKE_FAIL_STDERR" >&2
    exit "${FAKE_FAIL_CODE:-1}"
fi
exit 0
"""


@pytest.fixture
def fake_ersilia(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    script = bin_dir / "ersilia"
    script.write_text(FAKE_ERSILIA)
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    log = tmp_path / "calls.log"
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_ERSILIA_LOG", str(log))
    monkeypatch.delenv("FAKE_FAIL_CMD", raising=False)
    monkeypatch.delenv("FAKE_CLOSE_FAIL", raising=False)
    monkeypatch.chdir(tmp_path)

    def fail(cmd, stderr, returncode=1):
        monkeypatch.setenv("FAKE_FAIL_CMD", cmd)
        monkeypatch.setenv("FAKE_FAIL_STDERR", stderr)
        monkeypatch.setenv("FAKE_FAIL_CODE", str(returncode))

    def calls():
        return log.read_text().split() if log.exists() else []

    return fail, calls


@pytest.fixture
def inspector(tmp_path, monkeypatch):
    monkeypatch.setattr(
        inspect_module.IOService, "get_model_type", staticmethod(lambda *a: None)
    )
    return ModelInspector("eos0abc", str(tmp_path))


def test_successful_run_records_time_and_closes(fake_ersilia, inspector):
    _, calls = fake_ersilia
    result, timeout = inspector._run_performance_check(1, False)
    assert result.success
    assert not timeout
    assert re.fullmatch(PERFORMANCE_PATTERN, result.details)
    assert calls() == ["serve", "example", "run", "close"]


@pytest.mark.parametrize(
    "message",
    [
        "The input could not be processed.",
        "None of the 10,000 inputs could be processed.",
    ],
)
def test_no_input_processed_is_not_fatal(fake_ersilia, inspector, message):
    fail, calls = fake_ersilia
    fail("run", f"  ✖  {message}")
    result, timeout = inspector._run_performance_check(1, False)
    assert result.success
    assert not timeout
    assert re.fullmatch(PERFORMANCE_PATTERN, result.details)
    assert calls() == ["serve", "example", "run", "close"]


@pytest.mark.parametrize("step", ["serve", "example", "run", "close"])
def test_other_failures_still_fail_and_close(fake_ersilia, inspector, step):
    fail, calls = fake_ersilia
    fail(step, "  ✖  Model eos0abc could not be served.")
    result, timeout = inspector._run_performance_check(1, False)
    assert not result.success
    assert not timeout
    assert result.details == "Something happened when running the deep check!"
    assert calls()[-1] == "close"


@pytest.mark.parametrize("step", ["serve", "example", "close"])
def test_empty_input_message_outside_run_still_fails(fake_ersilia, inspector, step):
    fail, calls = fake_ersilia
    fail(step, "The input could not be processed.")
    result, timeout = inspector._run_performance_check(1, False)
    assert not result.success
    assert not timeout
    assert result.details == "Something happened when running the deep check!"
    assert calls()[-1] == "close"


def test_close_failure_after_empty_run_still_fails(fake_ersilia, inspector, monkeypatch):
    fail, calls = fake_ersilia
    fail("run", "The input could not be processed.")
    monkeypatch.setenv("FAKE_CLOSE_FAIL", "1")
    result, timeout = inspector._run_performance_check(1, False)
    assert not result.success
    assert not timeout
    assert result.details == "Something happened when running the deep check!"
    assert calls() == ["serve", "example", "run", "close"]


def test_other_exit_code_with_empty_input_message_still_fails(fake_ersilia, inspector):
    fail, calls = fake_ersilia
    fail("run", "The input could not be processed.", returncode=2)
    result, timeout = inspector._run_performance_check(1, False)
    assert not result.success
    assert not timeout
    assert result.details == "Something happened when running the deep check!"
    assert calls() == ["serve", "example", "run", "close"]


def test_deep_check_passes_when_example_inputs_return_nothing(fake_ersilia, inspector):
    fail, calls = fake_ersilia
    fail("run", "  ✖  The input could not be processed.")
    result = inspector.check_computational_performance()
    assert result.success
    timings = re.findall(PERFORMANCE_PATTERN, result.details)
    assert [int(n) for n, _ in timings] == [1, 10, 100, 1000, 10000]
    assert calls().count("serve") == calls().count("close") == 5
