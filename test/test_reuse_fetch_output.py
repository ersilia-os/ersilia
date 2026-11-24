"""Tests for reusing the output that 'ersilia fetch' produced for the example input."""

from unittest.mock import MagicMock

import pytest

from ersilia.default import (
    EXAMPLE_STANDARD_INPUT_CSV_FILENAME,
    EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME,
    PREDEFINED_COLUMN_FILE,
    PREDEFINED_EXAMPLE_INPUT_FILES,
)
from ersilia.publish.test.services.checks import CheckService
from ersilia.publish.test.services.constants import STATUS_CONFIGS, Options
from ersilia.publish.test.services.runner import RunnerService, _fetched_example_output

PASSED = str(STATUS_CONFIGS.PASSED)
FAILED = str(STATUS_CONFIGS.FAILED)

EXAMPLE_INPUT = "smiles\nCCO\nCCN\nc1ccccc1\n"
EXAMPLE_OUTPUT = "key,input,score\nk1,CCO,0.1\nk2,CCN,0.2\nk3,c1ccccc1,0.3\n"


@pytest.fixture
def dirs(tmp_path):
    model_path = tmp_path / "dest"
    repo = tmp_path / "repo"
    model_path.mkdir()
    repo.mkdir()
    (model_path / EXAMPLE_STANDARD_INPUT_CSV_FILENAME).write_text(EXAMPLE_INPUT)
    (model_path / EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME).write_text(EXAMPLE_OUTPUT)
    test_input = repo / "run_input.csv"
    test_input.write_text(EXAMPLE_INPUT)
    return model_path, test_input


def test_reuses_output_when_inputs_match(dirs):
    model_path, test_input = dirs
    assert _fetched_example_output(str(model_path), str(test_input)) == str(
        model_path / EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME
    )


def test_no_reuse_when_inputs_differ(dirs):
    # The image's run_input.csv can be older than the repository's.
    model_path, test_input = dirs
    test_input.write_text("smiles\nCCO\nCCN\nCCCl\n")
    assert _fetched_example_output(str(model_path), str(test_input)) is None


@pytest.mark.parametrize(
    "missing", [EXAMPLE_STANDARD_INPUT_CSV_FILENAME, EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME]
)
def test_no_reuse_when_fetch_files_are_missing(dirs, missing):
    model_path, test_input = dirs
    (model_path / missing).unlink()
    assert _fetched_example_output(str(model_path), str(test_input)) is None


def test_no_reuse_when_fetch_output_is_empty(dirs):
    model_path, test_input = dirs
    (model_path / EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME).write_text("")
    assert _fetched_example_output(str(model_path), str(test_input)) is None


@pytest.mark.parametrize("test_input", [None, "does/not/exist.csv"])
def test_no_reuse_without_a_test_input(dirs, test_input):
    model_path, _ = dirs
    assert _fetched_example_output(str(model_path), test_input) is None


def _check_service(repo):
    service = CheckService.__new__(CheckService)
    service.logger = MagicMock()
    service.dir = str(repo)
    service.ios = MagicMock()
    service.ios.get_output_consistency.return_value = "Variable"
    return service


@pytest.fixture
def repo(tmp_path):
    repo = tmp_path / "repo"
    (repo / PREDEFINED_EXAMPLE_INPUT_FILES[0]).parent.mkdir(parents=True)
    (repo / PREDEFINED_COLUMN_FILE).parent.mkdir(parents=True)
    (repo / PREDEFINED_EXAMPLE_INPUT_FILES[0]).write_text(EXAMPLE_INPUT)
    (repo / PREDEFINED_COLUMN_FILE).write_text("name,type\nscore,float\n")
    return repo


def test_simple_run_checks_the_fetch_output_without_running(repo, tmp_path):
    fetch_output = tmp_path / "fetch_output.csv"
    fetch_output.write_text(EXAMPLE_OUTPUT)
    run_model = MagicMock()

    status = _check_service(repo).check_simple_model_output(
        run_model, reuse_output=str(fetch_output)
    )

    run_model.assert_not_called()
    with open(Options.OUTPUT_CSV.value) as f:
        assert f.read() == EXAMPLE_OUTPUT
    assert [s[2] for s in status] == [PASSED, PASSED]


def test_simple_run_still_fails_a_bad_fetch_output(repo, tmp_path):
    # Reusing the output must not hide a problem the check would catch.
    fetch_output = tmp_path / "fetch_output.csv"
    fetch_output.write_text("key,input,score\nk1,CCO,\nk2,CCN,\nk3,c1ccccc1,\n")

    status = _check_service(repo).check_simple_model_output(
        MagicMock(), reuse_output=str(fetch_output)
    )

    assert status[0][2] == FAILED


def test_simple_run_runs_the_model_without_a_fetch_output(repo):
    def run_model(inputs, output, batch):
        with open(output, "w") as f:
            f.write(EXAMPLE_OUTPUT)

    run_model = MagicMock(side_effect=run_model)

    status = _check_service(repo).check_simple_model_output(run_model)

    run_model.assert_called_once()
    assert run_model.call_args.kwargs["inputs"] == str(
        repo / PREDEFINED_EXAMPLE_INPUT_FILES[0]
    )
    assert [s[2] for s in status] == [PASSED, PASSED]


def test_runner_finds_the_fetch_output_for_the_repository_input(repo, tmp_path):
    model_path = tmp_path / "dest"
    model_path.mkdir()
    (model_path / EXAMPLE_STANDARD_INPUT_CSV_FILENAME).write_text(EXAMPLE_INPUT)
    (model_path / EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME).write_text(EXAMPLE_OUTPUT)
    runner = RunnerService.__new__(RunnerService)
    runner.logger = MagicMock()
    runner.dir = str(repo)
    runner.model_path = str(model_path)

    assert runner._reusable_fetch_output() == str(
        model_path / EXAMPLE_STANDARD_OUTPUT_CSV_FILENAME
    )
    runner.dir = None
    assert runner._reusable_fetch_output() is None
