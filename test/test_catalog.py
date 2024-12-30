import json
import os
import pytest
from ersilia.hub.content.catalog import CatalogTable


@pytest.fixture
def catalog_samples():
    file_path = os.path.join(
        os.path.dirname(__file__), "inputs", "catalog_samples.json"
    )
    with open(file_path, "r") as f:
        samples = json.load(f)
    return samples


def test_as_list_of_dicts(catalog_samples):
    columns = ["Identifier", "Slug", "Title"]

    # Test with standard catalog samples
    catalog_table = CatalogTable(
        data=[list(item.values()) for item in catalog_samples], columns=columns
    )
    result = catalog_table.as_list_of_dicts()
    assert (
        result == catalog_samples
    ), "The result does not match the expected catalog samples"

    # Test with empty catalog data
    catalog_table_empty = CatalogTable(data=[], columns=columns)
    result_empty = catalog_table_empty.as_list_of_dicts()
    assert result_empty == [], "The result should be an empty list for empty input data"
