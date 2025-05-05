import pytest, asyncio
from ersilia.core.model import ErsiliaModel
from ersilia.hub.fetch.fetch import ModelFetcher

from ..cli.utils import create_compound_input_csv
import csv

INPUT_CSV = "inputs.csv"
OUTPUT_CSV = "outputs.csv"
MODEL_ID = "eos3b5e"
MODEL_EXIST_MESSAGE = "Model already exists on your system. If you want to fetch it again, please delete the existing model first."

def fetch():
    mf = ModelFetcher(
            overwrite=True,
            force_from_github=True,
        )
    res = asyncio.run(mf.fetch(MODEL_ID))
    return res

def simple_csv_content_check(output):
    checks = []
    with open(output, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                checks.append(False)
    checks = all(chk for chk in checks)
    return checks


@pytest.fixture
def compound_csv():
    create_compound_input_csv(INPUT_CSV)
    yield INPUT_CSV

def test_api(compound_csv):
    is_fetched = fetch()
    model = ErsiliaModel(model=MODEL_ID, verbose=True, output_source=None)
    model.serve()
    model.run(input=INPUT_CSV, output=OUTPUT_CSV)
    has_contents = simple_csv_content_check(OUTPUT_CSV)

    assert is_fetched.fetch_success is True or MODEL_EXIST_MESSAGE in is_fetched.reason
    assert has_contents is True