import json
import os
import pytest
from ersilia.hub.content.catalog import CatalogTable


@pytest.fixture
def catalog_samples():
    """
    Fixture that loads test data from a JSON file.
    The JSON file contains realtime catalog samples for testing
    """
    file_path = os.path.join(os.path.dirname(__file__), 'inputs', 'test_catalog','catalog_sample.json')
    with open(file_path, 'r') as f:
        samples = json.load(f)
    return samples


@pytest.fixture
def catalog_columns():
    """Fixture that provides coloumn names used in testing"""
    return ['Identifier', 'Slug', 'Title']