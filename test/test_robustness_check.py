import logging
from collections import namedtuple

import pytest

from ersilia.default import NO_INPUT_PROCESSED_EXIT_CODE
from ersilia.publish.test.services import constants
from ersilia.publish.test.services.checks import CheckService
from ersilia.publish.test.services.constants import STATUS_CONFIGS, Options

CommandResult = namedtuple("CommandResult", ["returncode", "stdout", "stderr"])

# One planned outcome per attempt (3, 5, then 10 samples)
NO_INPUT = "no input processed"  # exit code 3, no output file
CRASH = "crash"  # exit code 1, no output file
RAISES = "raises"  # run_model raises, no output file
NULL = "null"  # output file with only empty values
VALUE = "value"  # output file with a value


def _service(output_consistency):
    service = CheckService.__new__(CheckService)
    service.logger = logging.getLogger("test_robustness_check")
    service._get_output_consistency = lambda: output_consistency
    return service


def _run(output_consistency, plan):
    attempts = iter(plan)
    calls = []

    def run_example(n_samples, file_name, simple, mode):
        calls.append(n_samples)

    def run_model(inputs, output, batch):
        outcome = next(attempts)
        if outcome == NO_INPUT:
            return CommandResult(NO_INPUT_PROCESSED_EXIT_CODE, "", "")
        if outcome == CRASH:
            return CommandResult(1, "", "")
        if outcome == RAISES:
            raise RuntimeError("serve failed")
        value = "CCO" if outcome == VALUE else ""
        with open(output, "w") as f:
            f.write(f"key,input,smi_0\nk,C,{value}\n")
        return CommandResult(0, "", "")

    result = _service(output_consistency).check_robustness(run_example, run_model)
    _, detail, status = result[0]
    return detail, status, calls


@pytest.fixture(autouse=True)
def _isolated_eos_tmp(tmp_path, monkeypatch):
    # Options.*.value resolves file names under EOS_TMP/files, which is shared
    # with real `ersilia test` runs; point it at a fresh folder per test
    monkeypatch.setattr(constants, "EOS_TMP", str(tmp_path))
    yield
    assert Options.ROBUSTNESS_OUTPUT_CSV.value.startswith(str(tmp_path))


PASSED = str(STATUS_CONFIGS.PASSED)
FAILED = str(STATUS_CONFIGS.FAILED)


def test_variable_model_moves_on_when_no_input_is_processed():
    detail, status, calls = _run("Variable", [NO_INPUT, NO_INPUT, VALUE])
    assert status == PASSED
    assert "10 sample" in detail
    assert calls == [3, 5, 10]


def test_variable_model_fails_when_no_input_is_ever_processed():
    detail, status, calls = _run("Variable", [NO_INPUT, NO_INPUT, NO_INPUT])
    assert status == FAILED
    assert "no input processed" in detail
    assert calls == [3, 5, 10]


def test_variable_model_fails_on_mixed_null_and_no_input():
    _, status, calls = _run("Variable", [NULL, NO_INPUT, NULL])
    assert status == FAILED
    assert calls == [3, 5, 10]


@pytest.mark.parametrize("outcome", [CRASH, RAISES])
def test_variable_model_fails_at_once_on_a_crash(outcome):
    detail, status, calls = _run("Variable", [outcome, VALUE, VALUE])
    assert status == FAILED
    assert "Output file not generated with 3 samples" in detail
    assert calls == [3]


def test_stale_output_file_is_not_read_as_a_result():
    with open(Options.ROBUSTNESS_OUTPUT_CSV.value, "w") as f:
        f.write("key,input,smi_0\nk,C,CCO\n")
    _, status, _ = _run("Variable", [NO_INPUT, NO_INPUT, NO_INPUT])
    assert status == FAILED


def test_fixed_model_fails_at_once_when_no_input_is_processed():
    detail, status, calls = _run("Fixed", [NO_INPUT])
    assert status == FAILED
    assert "exit code 3" in detail
    assert calls == [3]


def test_fixed_model_passes_with_a_value():
    _, status, calls = _run("Fixed", [VALUE])
    assert status == PASSED
    assert calls == [3]
