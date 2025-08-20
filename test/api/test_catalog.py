from unittest.mock import patch, ANY
import pytest
import pandas as pd

from ersilia.api.create_api import ErsiliaCatalog
from ersilia.api.commands.catalog import catalog


@pytest.fixture
def mock_catalog():
    # Patch the function accessed as `ersilia.api.create_api.catalog.catalog(...)`
    with patch("ersilia.api.create_api.catalog.catalog") as mock_catalog_func:
        mock_df = pd.DataFrame({
            "Identifier": ["eos3b5e", "eos78ao"],
            "Slug": ["molecular-weight", "mordred"]
        })
        mock_catalog_func.return_value = mock_df
        yield mock_catalog_func


def test_ersilia_catalog_hub_true(mock_catalog):
    catalog_obj = ErsiliaCatalog()
    result = catalog_obj.catalog(hub=True)
    mock_catalog.assert_called_once_with(
        hub=True,
        file_name=ANY, browser=ANY, more=ANY, card=ANY,
        model=ANY, as_json=ANY, verbose=ANY
    )
    assert result is not None


def test_ersilia_catalog_hub_false(mock_catalog):
    catalog_obj = ErsiliaCatalog()
    result = catalog_obj.catalog(hub=False)
    mock_catalog.assert_called_once_with(
        hub=False,
        file_name=ANY, browser=ANY, more=ANY, card=ANY,
        model=ANY, as_json=ANY, verbose=ANY
    )
    assert result is not None


def test_ersilia_catalog_default(mock_catalog):
    catalog_obj = ErsiliaCatalog()
    result = catalog_obj.catalog()
    mock_catalog.assert_called_once_with(
        hub=False,  # default
        file_name=ANY, browser=ANY, more=ANY, card=ANY,
        model=ANY, as_json=ANY, verbose=ANY
    )
    assert result is not None


def test_catalog_function_hub_parameter():
    import inspect
    sig = inspect.signature(catalog)
    assert 'hub' in sig.parameters
    assert sig.parameters['hub'].default is False
