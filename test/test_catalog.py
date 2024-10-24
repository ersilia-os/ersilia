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

def test_as_table(catalog_samples, catalog_columns):
    """
    Test the as_table method of CatalogTable

    Verifes that the method produces correctly formatted table string
    with proper borders and content alignment
    """

    data = []
    for sample in catalog_samples:
        row = [sample.get(col) for col in catalog_columns]
        data.append(row)
    
    catalog_table = CatalogTable(data=data, columns=catalog_columns)
    table_result = catalog_table.as_table()
    
    # Basic Structure Test
    assert isinstance(table_result, str), "Table result should be a string"
    lines = table_result.split('\n')
    assert len(lines) > 0, "Table should not be empty"
    
    # Check Table Structure
    assert lines[0].startswith('┌'), "Table should start with top-left corner"
    assert lines[0].endswith('┐'), "Table should end with top-right corner"
    
    # Check the HeaderRow
    header_row = lines[1]
    assert header_row.startswith('│'), "Header should start with vertical border"
    assert header_row.endswith('│'), "Header should end with vertical border"
    for column in catalog_columns:
        assert column in header_row, f"Header should contain column '{column}'"
    
    # Check data rows
    data_start_index = 3  
    for i, sample in enumerate(catalog_samples):
        row_index = data_start_index + (i * 2)  # Account for lines in between data '|'
        if row_index < len(lines):
            data_row = lines[row_index]
            assert data_row.startswith('│'), f"Data row {i} should start with vertical border"
            assert data_row.endswith('│'), f"Data row {i} should end with vertical border"
            # Check if sample data is present in the row
            assert sample['Identifier'] in data_row, f"Row {i} should contain Identifier"
            assert sample['Slug'] in data_row, f"Row {i} should contain Slug"
    
    # Check bottom 
    assert lines[-1].startswith('└'), "Table should end with bottom-left corner"
    assert lines[-1].endswith('┘'), "Table should end with bottom-right corner"
    
    # Test empty table
    catalog_table_empty = CatalogTable(data=[], columns=catalog_columns)
    empty_table_result = catalog_table_empty.as_table()
    empty_lines = empty_table_result.split('\n')
    
    assert len(empty_lines) == 4, "Empty table should have 4 lines (top border, header, separator, bottom border)"
    assert empty_lines[0].startswith('┌'), "Empty table should start with top-left corner"
    assert empty_lines[-1].startswith('└'), "Empty table should end with bottom-left corner"
    
    # Test table with None values
    data_with_none = [[None if j == 1 else val for j, val in enumerate(row)] for row in data]
    catalog_table_none = CatalogTable(data=data_with_none, columns=catalog_columns)
    none_table_result = catalog_table_none.as_table()
    
    assert isinstance(none_table_result, str), "Table with None values should be a string"
    none_lines = none_table_result.split('\n')
    assert len(none_lines) > 0, "Table with None values should not be empty"
    
    # Test column width
    max_width = max(
        max(len(str(item)) if item is not None else 0 for item in [col] + [row[i] for row in data])
        for i, col in enumerate(catalog_columns)
    )
    assert any(len(line) >= max_width for line in lines), "Table should accommodate longest content"
