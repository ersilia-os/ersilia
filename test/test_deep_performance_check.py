import subprocess
from unittest.mock import patch

from ersilia.default import NO_INPUT_PROCESSED_EXIT_CODE
from ersilia.publish.test.services.inspect import ModelInspector


def _inspector():
    inspector = ModelInspector.__new__(ModelInspector)
    inspector.model = "eos4qda"
    return inspector


def _completed(returncode):
    return subprocess.CompletedProcess(args="", returncode=returncode)


def test_inputs_that_could_not_be_processed_still_record_the_time():
    with patch("subprocess.run", return_value=_completed(NO_INPUT_PROCESSED_EXIT_CODE)) as run:
        result, timeout = _inspector()._run_performance_check(1, False)
    assert result.success
    assert "1 predictions executed in" in result.details
    assert not timeout
    # The chain stopped at 'ersilia run', so the model is closed separately.
    assert run.call_args_list[-1].args[0] == "ersilia close"


def test_other_failures_still_fail_the_check():
    with patch("subprocess.run", return_value=_completed(1)):
        result, timeout = _inspector()._run_performance_check(1, False)
    assert not result.success
    assert not timeout
