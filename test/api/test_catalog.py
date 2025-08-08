from unittest.mock import MagicMock, patch
import pytest
import pandas as pd

from ersilia.api.create_api import ErsiliaCatalog
from ersilia.api.commands.catalog import catalog


# Fixture to mock the catalog function
@pytest.fixture
def mock_catalog():
    with patch("ersilia.api.commands.catalog.catalog") as mock_catalog_func:
        # Create a mock DataFrame for testing
        mock_df = pd.DataFrame({
            "Identifier": ["eos3b5e", "eos78ao"],
            "Slug": ["molecular-weight", "mordred"]
        })
        mock_catalog_func.return_value = mock_df
        yield mock_catalog_func


# Test ErsiliaCatalog.catalog() with hub=True
def test_ersilia_catalog_hub_true(mock_catalog):
    """Test that ErsiliaCatalog.catalog(hub=True) calls the underlying catalog function with hub=True."""
    catalog_obj = ErsiliaCatalog()
    result = catalog_obj.catalog(hub=True)
    
    # Verify the underlying catalog function was called with hub=True
    mock_catalog.assert_called_once_with(hub=True)
    assert result is not None


# Test ErsiliaCatalog.catalog() with hub=False (default)
def test_ersilia_catalog_hub_false(mock_catalog):
    """Test that ErsiliaCatalog.catalog(hub=False) calls the underlying catalog function with hub=False."""
    catalog_obj = ErsiliaCatalog()
    result = catalog_obj.catalog(hub=False)
    
    # Verify the underlying catalog function was called with hub=False
    mock_catalog.assert_called_once_with(hub=False)
    assert result is not None


# Test ErsiliaCatalog.catalog() with default parameter
def test_ersilia_catalog_default(mock_catalog):
    """Test that ErsiliaCatalog.catalog() calls the underlying catalog function with default hub=False."""
    catalog_obj = ErsiliaCatalog()
    result = catalog_obj.catalog()
    
    # Verify the underlying catalog function was called with default hub=False
    mock_catalog.assert_called_once_with(hub=False)
    assert result is not None


# Test that the catalog function accepts the hub parameter
def test_catalog_function_hub_parameter():
    """Test that the catalog function accepts the hub parameter."""
    # This test verifies that the catalog function signature includes the hub parameter
    import inspect
    sig = inspect.signature(catalog)
    assert 'hub' in sig.parameters
    assert sig.parameters['hub'].default is False 