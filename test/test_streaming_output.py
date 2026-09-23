import h5py
import numpy as np
import pytest

from ersilia.io.output import GenericOutputAdapter

N_ROWS = 25
CHUNK = 10


def _results(dtype):
    rows = []
    for i in range(N_ROWS):
        if dtype == "float":
            out = {
                f"f{j}": (None if (i + j) % 7 == 0 else i * 0.5 + j) for j in range(4)
            }
        else:
            out = {f"f{j}": f"s{i}_{j}" for j in range(4)}
        rows.append({"input": {"key": f"k{i}", "input": f"C{i}"}, "output": out})
    return rows


def _adapter(dtype):
    columns_info = {"name": [f"f{j}" for j in range(4)], "type": [dtype] * 4}
    return GenericOutputAdapter(model_id="eos0xxx", columns_info=columns_info)


def _write_whole(adapter, results, path):
    df = adapter._to_dataframe(results)
    df.write(str(path), delimiter="\t" if str(path).endswith("tsv") else ",")


def _write_chunked(adapter, results, path):
    for start in range(0, len(results), CHUNK):
        adapter.write_chunk(results[start : start + CHUNK], str(path), start > 0)


@pytest.mark.parametrize("dtype", ["float", "string"])
@pytest.mark.parametrize("ext", ["csv", "tsv"])
def test_chunked_text_matches_single_write(tmp_path, dtype, ext):
    """Streaming batches to CSV/TSV gives the same file as one full write."""
    adapter = _adapter(dtype)
    results = _results(dtype)
    _write_whole(adapter, results, tmp_path / f"whole.{ext}")
    _write_chunked(adapter, results, tmp_path / f"chunked.{ext}")
    whole = (tmp_path / f"whole.{ext}").read_text()
    assert whole == (tmp_path / f"chunked.{ext}").read_text()
    assert len(whole.splitlines()) == N_ROWS + 1


@pytest.mark.parametrize("dtype", ["float", "string"])
def test_chunked_hdf5_matches_single_write(tmp_path, dtype):
    """Streaming batches to HDF5 gives the same datasets as one full write."""
    adapter = _adapter(dtype)
    results = _results(dtype)
    _write_whole(adapter, results, tmp_path / "whole.h5")
    _write_chunked(adapter, results, tmp_path / "chunked.h5")
    with h5py.File(tmp_path / "whole.h5") as a, h5py.File(tmp_path / "chunked.h5") as b:
        assert (
            set(a.keys()) == set(b.keys()) == {"Values", "Keys", "Inputs", "Features"}
        )
        for name in a:
            x, y = a[name][:], b[name][:]
            assert x.shape == y.shape and x.dtype == y.dtype
            if x.dtype.kind == "f":
                assert np.array_equal(x, y, equal_nan=True)
            else:
                assert np.array_equal(x, y)
        assert b["Values"].shape == (N_ROWS, 4)
        assert b["Features"].shape == (4,)
