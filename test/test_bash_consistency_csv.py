"""Tests for the CSV reading and comparison behind the bash/Ersilia consistency check."""

import pytest

from ersilia.publish.test.services.constants import STATUS_CONFIGS
from ersilia.publish.test.services.runner import (
    _compare_outputs,
    _compute_rmse,
    _read_output_csv,
    _write_ersilia_input,
)

PASSED = str(STATUS_CONFIGS.PASSED)
FAILED = str(STATUS_CONFIGS.FAILED)

CXSMILES = "COC(C)(CNC(=O)[C@H]1C[C@@H]1C1=CC=CC=C1Br)C1CCCCC1 |&1:8,10|"


def _exact(a, b):
    return a == b


def _write(path, text):
    path.write_text(text)
    return str(path)


def test_quoted_cell_with_commas_stays_in_one_column(tmp_path):
    path = _write(
        tmp_path / "out.csv",
        f'smiles_00,score_00,smiles_01\n"{CXSMILES}",0.91,CCO\n',
    )
    rows, error = _read_output_csv(path)
    assert error is None
    assert rows == [{"smiles_00": CXSMILES, "score_00": 0.91, "smiles_01": "CCO"}]


def test_skip_columns_drops_key_and_input(tmp_path):
    path = _write(
        tmp_path / "out.csv",
        f'key,input,smiles_00,score_00\nabc,"C,C","{CXSMILES}",7\n',
    )
    rows, error = _read_output_csv(path, skip_columns=2)
    assert error is None
    assert rows == [{"smiles_00": CXSMILES, "score_00": 7}]


def test_empty_file_reports_error(tmp_path):
    rows, error = _read_output_csv(_write(tmp_path / "out.csv", ""))
    assert rows == []
    assert error == "File is empty"


def test_cxsmiles_outputs_are_compared_not_skipped(tmp_path):
    bsh, _ = _read_output_csv(
        _write(
            tmp_path / "bash.csv",
            f'smiles_00,score_00\n"{CXSMILES}",0.5\nCCO,0.4\n',
        )
    )
    ers, _ = _read_output_csv(
        _write(
            tmp_path / "ersilia.csv",
            f'key,input,smiles_00,score_00\nk1,CC,"{CXSMILES}",0.5\nk2,CO,CCO,0.4\n',
        ),
        skip_columns=2,
    )
    status = _compare_outputs(bsh, ers, _exact)
    assert ("String Similarity", "> 95%", PASSED) in status
    assert all(s[2] == PASSED for s in status)


def test_differing_cxsmiles_fails_string_similarity():
    bsh = [{"smiles_00": CXSMILES}]
    ers = [{"smiles_00": "CCO |&1:8,10|"}]
    status = _compare_outputs(bsh, ers, _exact)
    assert status == [("String Similarity", "< 95%", FAILED)]


def test_mixed_type_column_fails():
    bsh = [{"col": "CCO"}, {"col": 12}]
    ers = [{"col": "CCO"}, {"col": 12}]
    status = _compare_outputs(bsh, ers, _exact)
    assert len(status) == 1
    name, detail, outcome = status[0]
    assert name == "Type Check-col"
    assert "mixes numbers and text" in detail
    assert outcome == FAILED


def test_numeric_column_with_blanks_passes():
    bsh = [{"score": 0.5}, {"score": ""}, {"score": 1}]
    ers = [{"score": 0.5}, {"score": ""}, {"score": 1.0}]
    status = _compare_outputs(bsh, ers, _exact)
    assert status == [("RMSE-MEAN", "RMSE < 10% | 0.0%", PASSED)]


def test_compute_rmse_divides_inside_square_root():
    assert _compute_rmse([0.0] * 10, [0.3] * 10) == pytest.approx(0.3)
    assert _compute_rmse([0.0, 3.0], [4.0, 3.0]) == pytest.approx((16 / 2) ** 0.5)


def test_inconsistent_numeric_column_fails_rmse_for_many_rows():
    bsh = [{"score": 0.0} for _ in range(100)]
    ers = [{"score": 1.0} for _ in range(100)]
    status = _compare_outputs(bsh, ers, _exact)
    assert status == [("RMSE-score", "RMSE > 10%: 100.0%", FAILED)]


def test_string_column_with_blanks_passes():
    bsh = [{"smiles": "CCO"}, {"smiles": ""}]
    ers = [{"smiles": "CCO"}, {"smiles": ""}]
    status = _compare_outputs(bsh, ers, _exact)
    assert ("String Similarity", "> 95%", PASSED) in status
    assert all(s[2] == PASSED for s in status)


def test_type_mismatch_message_shows_both_values():
    bsh = [{"score": "abc"}]
    ers = [{"score": 0.5}]
    status = _compare_outputs(bsh, ers, _exact)
    name, detail, outcome = status[0]
    assert outcome == FAILED
    assert "bash value='abc' (str)" in detail
    assert "ersilia value=0.5 (float)" in detail


@pytest.mark.parametrize(
    "bsh, ers",
    [([], [{"a": 1}]), ([{"a": 1}], []), ([], [])],
)
def test_empty_output_fails(bsh, ers):
    status = _compare_outputs(bsh, ers, _exact)
    assert len(status) == 1
    assert status[0][2] == FAILED
    assert "Cannot compare outputs" in status[0][1]


def test_no_shared_columns_fails():
    status = _compare_outputs([{"a": 1}], [{"b": 1}], _exact)
    assert status[0][2] == FAILED
    assert "share no column names" in status[0][1]


def test_bash_input_mirrors_what_ersilia_passes_to_run_sh(tmp_path):
    example = tmp_path / "run_input.csv"
    example.write_bytes(f'smiles\nCCO\n"{CXSMILES}"\n\nc1ccccc1\n'.encode())
    path = tmp_path / "example_file.csv"

    _write_ersilia_input(str(example), str(path))

    # ersilia-pack writes an 'input' header, one input per row, via csv.writer.
    assert path.read_bytes() == (
        f'input\r\nCCO\r\n"{CXSMILES}"\r\nc1ccccc1\r\n'.encode()
    )


@pytest.mark.parametrize("line_end", ["\n", "\r\n"])
def test_bash_input_leaves_the_model_example_untouched(tmp_path, line_end):
    example = tmp_path / "run_input.csv"
    original = line_end.join(["input", "CCO", "CCN", ""]).encode()
    example.write_bytes(original)

    _write_ersilia_input(str(example), str(tmp_path / "example_file.csv"))

    assert example.read_bytes() == original
