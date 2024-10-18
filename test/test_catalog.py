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
    file_path = os.path.join(os.path.dirname(__file__), 'inputs', 'test_catalog', 'catalog_samples.json')
    
    assert os.path.exists(file_path), f"File not found: {file_path}"
    
    with open(file_path, 'r') as f:
        data = json.load(f)
    samples = data.get('catalog_samples', [])
    assert isinstance(samples, list), "Expected a list of catalog samples"
    assert all(isinstance(item, dict) for item in samples), "Each sample should be a dictionary"
    return samples

@pytest.fixture
def catalog_columns():
    """
    Fixture that provides column names used in testing.
    Returns the columns we want to test with.
    """
    return ['Identifier', 'Slug', 'Title', 'Description', 'Status', 'Task', 
            'Input', 'Output', 'Tag', 'Publication', 'License', 'Year', 
            'Dependencies', 'Contributors']

