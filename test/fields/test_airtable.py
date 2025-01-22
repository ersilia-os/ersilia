import sys
import os
import pytest
from unittest.mock import patch, MagicMock

MODEL_ID = "eos3b5e"
ENVIRONMENT_SIZE = 345
IMAGE_SIZE = 1010
COMPUTATIONAL_PERFORMANCE_1 = 10.11
COMPUTATIONAL_PERFORMANCE_10 = 14.12
COMPUTATIONAL_PERFORMANCE_100 = 50.45
DUMMP_API_KEY = "airtable_api_key_456"

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.github/scripts")))

from airtableops import AirtableInterface, AirtableMetadata

MOCK_RECORD = {
    "fields": {
        "Identifier": MODEL_ID,
        "Environment Size": ENVIRONMENT_SIZE,
        "Image Size": IMAGE_SIZE,
        "Computational Performance 1": COMPUTATIONAL_PERFORMANCE_1,
        "Computational Performance 10": COMPUTATIONAL_PERFORMANCE_10,
        "Computational Performance 100": COMPUTATIONAL_PERFORMANCE_100,
    },
    "id": "rec_123",
}


@pytest.fixture
def mock_ro_api_key():
    return DUMMP_API_KEY


@pytest.fixture
def mock_rw_api_key():
    return DUMMP_API_KEY


@pytest.fixture
def mock_airtable_table():
    mock_table = MagicMock()
    mock_table.iterate.return_value = [[MOCK_RECORD]]
    mock_table.all.return_value = [MOCK_RECORD]
    return mock_table


@patch("airtableops.AirtableInterface._get_ro_airtable_api_key")
@patch("airtableops.AirtableInterface._create_table")
def test_airtable_interface_ro(mock_create_table, mock_get_ro_api_key, mock_ro_api_key, mock_airtable_table):
    mock_get_ro_api_key.return_value = mock_ro_api_key
    mock_create_table.return_value = mock_airtable_table

    ai = AirtableInterface(mode="ro")
    items = list(ai.items())
    assert len(items) == 1
    assert items[0]["fields"]["Identifier"] == MODEL_ID


@patch("airtableops.AirtableInterface._create_table")
def test_airtable_interface_rw(mock_create_table, mock_rw_api_key, mock_airtable_table):
    mock_create_table.return_value = mock_airtable_table

    ai = AirtableInterface(mode="rw", api_key=mock_rw_api_key)
    items = list(ai.items_all())
    assert len(items) == 1
    assert items[0]["fields"]["Identifier"] == MODEL_ID
    assert items[0]["fields"]["Environment Size"] == ENVIRONMENT_SIZE
    assert items[0]["fields"]["Image Size"] == IMAGE_SIZE
    assert items[0]["fields"]["Computational Performance 1"] == COMPUTATIONAL_PERFORMANCE_1
    assert items[0]["fields"]["Computational Performance 10"] == COMPUTATIONAL_PERFORMANCE_10
    assert items[0]["fields"]["Computational Performance 100"] == COMPUTATIONAL_PERFORMANCE_100


@patch("airtableops.AirtableMetadata._find_record")
def test_airtable_metadata_find_record(mock_find_record):
    mock_find_record.return_value = MOCK_RECORD["fields"]
    am = AirtableMetadata(model_id=MODEL_ID)
    data = am._find_record()
    assert data["Identifier"] == MODEL_ID
    assert data["Environment Size"] == ENVIRONMENT_SIZE
    assert data["Image Size"] == IMAGE_SIZE
    assert data["Computational Performance 1"] == COMPUTATIONAL_PERFORMANCE_1
    assert data["Computational Performance 10"] == COMPUTATIONAL_PERFORMANCE_10
    assert data["Computational Performance 100"] == COMPUTATIONAL_PERFORMANCE_100
