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

def test_as_json(catalog_samples, catalog_columns):
    """
    Test the as_json method of CatalogTable.
    
    Verifies that the method correctly converts table data to JSON format
    and maintains data integrity through the conversion.
    """

    #test standard samples
    data = []
    for sample in catalog_samples:
        row = [sample.get(col) for col in catalog_columns]
        data.append(row)
    
    catalog_table = CatalogTable(data=data, columns=catalog_columns)
    
    json_result = catalog_table.as_json()
    
    assert isinstance(json_result, str), "Result should be a string"
    
    parsed_result = json.loads(json_result)
    
    assert isinstance(parsed_result, list), "Parsed result should be a list"
    assert len(parsed_result) == len(catalog_samples), "Result should have same number of entries"
    
    for i, result_item in enumerate(parsed_result):
        assert isinstance(result_item, dict), "Each item should be a dictionary"
        assert set(result_item.keys()) == set(catalog_columns), "Keys should match columns"
        
        for col in catalog_columns:
            assert result_item[col] == catalog_samples[i][col], \
                f"Mismatch in {col} for item {i}"
    
    # Test empty catalog data
    catalog_table_empty = CatalogTable(data=[], columns=catalog_columns)
    json_result_empty = catalog_table_empty.as_json()
    parsed_result_empty = json.loads(json_result_empty)
    assert parsed_result_empty == [], "Empty table should produce empty JSON array"