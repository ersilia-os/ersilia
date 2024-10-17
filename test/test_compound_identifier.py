import pytest
from ersilia.utils.identifiers.compound import CompoundIdentifier
from unittest.mock import patch

@pytest.fixture
def compound_identifier():
    return CompoundIdentifier()

@pytest.mark.parametrize("header", ["smiles", "input", "SMILES", "INPUT"])
def test_is_input_header_positive(compound_identifier, header):
    """Test that valid input headers return True."""
    assert compound_identifier.is_input_header(header) is True

@pytest.mark.parametrize("header", ["key", "inchiKey", "KEY", "INCHIKEY"])
def test_is_key_header_positive(compound_identifier, header):
    """Test that valid key headers return True."""
    assert compound_identifier.is_key_header(header) is True

@patch('ersilia.utils.identifiers.compound.CompoundIdentifier._pubchem_smiles_to_inchikey')
def test_is_smiles_positive_chem_none(mock_pubchem, compound_identifier):
    compound_identifier.Chem = None
    mock_pubchem.return_value = "InChIKey"
    assert compound_identifier._is_smiles("valid_smiles")
