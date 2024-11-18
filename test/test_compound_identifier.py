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
@pytest.mark.parametrize("header", ["output", "invalid", "InChI", "inchikey"])
def test_is_input_header_negative(compound_identifier, header):
    """Test that invalid input headers return False."""
    assert compound_identifier.is_input_header(header) is False
    
@pytest.mark.parametrize("header", ["key", "inchiKey", "KEY", "INCHIKEY"])
def test_is_key_header_positive(compound_identifier, header):
    """Test that valid key headers return True."""
    assert compound_identifier.is_key_header(header) is True
    
@pytest.mark.parametrize("header", ["id","smiles","inchi","input", "some_header", "random", "header", ""])
def test_is_key_header_negative(compound_identifier, header):
    assert not compound_identifier.is_key_header(header)

@patch('ersilia.utils.identifiers.compound.CompoundIdentifier._pubchem_smiles_to_inchikey')
def test_is_smiles_positive_chem_none(mock_pubchem, compound_identifier):
    compound_identifier.Chem = None
    mock_pubchem.return_value = "InChIKey"
    
# Test with a valid SMILES input
    smiles_string = 'CCO' #Ethanol SMILES
    assert compound_identifier._is_smiles(smiles_string) is True

