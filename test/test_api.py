import pytest
from ersilia.core.model import ErsiliaModel
from .cli.utils import create_compound_input_csv
import csv

INPUT_CSV = "input.csv"
OUTPUT_CSV = "output.csv"

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
    model = ErsiliaModel(model="eos3b5e", verbose=True, service_class="conda")
    model.fetch()
    model.serve()
    model.run(input=INPUT_CSV, output=OUTPUT_CSV)
    has_contents = simple_csv_content_check(OUTPUT_CSV)
    assert has_contents is True