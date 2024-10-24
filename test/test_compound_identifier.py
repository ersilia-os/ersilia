import pytest
from ersilia.utils.identifiers.compound import CompoundIdentifier

@pytest.fixture
def compound_identifier():
    return CompoundIdentifier()

def test_is_smiles_with_chem():
    identifier = CompoundIdentifier()
    result = identifier._is_smiles("C(C(=O)O)c1ccccc1")
    assert result == True
    
def test_is_smiles_without_chem():
    identifier = CompoundIdentifier(local=False)
    result = identifier._is_smiles("C(C(=O)O)c1ccccc1")
    assert result == "WLJVXDMOQOGPHL-UHFFFAOYSA-N"

@pytest.mark.parametrize("header", ["smiles", "input", "SMILES", "INPUT"])
def test_is_input_header_positive(compound_identifier, header):
    """Test that valid input headers return True."""
    assert compound_identifier.is_input_header(header) is True
    
@pytest.mark.parametrize("header", ["key", "inchiKey", "KEY", "INCHIKEY"])
def test_is_key_header_positive(compound_identifier, header):
    """Test that valid key headers return True."""
    assert compound_identifier.is_key_header(header) is True

@pytest.mark.parametrize("header", ["id","smiles","inchi","input", "some_header", "random", "header", ""])
def test_is_key_header_negative(compound_identifier, header):
    assert not compound_identifier.is_key_header(header)
    
