"""Integration test of the Python API: needs Docker and network access."""

import pytest

from ersilia.api import Model
from ersilia.utils.exceptions_utils.api_exceptions import ModelNotServedError

MODEL_ID = "eos3b5e"
INPUTS = ["CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O", "CCO"]


def test_model_lifecycle(tmp_path):
    mdl = Model(MODEL_ID)
    assert isinstance(mdl.fetch(), bool)
    assert mdl.is_fetched()
    assert mdl.fetch() is False  # already fetched

    assert "card" in mdl.info()
    assert len(mdl.example(3)) == 3

    served = mdl.serve()
    assert served["model_id"] == MODEL_ID and served["url"].startswith("http")

    df = mdl.run(INPUTS)
    assert df.shape[0] == len(INPUTS)
    out = tmp_path / "out.csv"
    assert mdl.run(INPUTS, output=str(out)) == str(out) and out.exists()

    mdl.close()
    with pytest.raises(ModelNotServedError):
        mdl.run(INPUTS)

    mdl.delete()
    assert not mdl.is_fetched()
