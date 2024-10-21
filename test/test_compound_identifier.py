import pytest
from ersilia.utils.identifiers.compound import CompoundIdentifier

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

@patch('requests.get')
def test_nci_smiles_to_inchikey_positive(mock_get, compound_identifier):
    """Test _nci_smiles_to_inchikey with a mocked positive response."""
    mock_response = mock_get.return_value
    mock_response.status_code = 200
    mock_response.text = "InChIKey=BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    inchikey = compound_identifier._nci_smiles_to_inchikey("CCO")
    assert inchikey == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

@patch('requests.get')
def test_nci_smiles_to_inchikey_negative(mock_get, compound_identifier):
    """Test _nci_smiles_to_inchikey with a mocked negative response."""
    mock_response = mock_get.return_value
    mock_response.status_code = 404  

    inchikey = compound_identifier._nci_smiles_to_inchikey("invalid_smiles")
    assert inchikey is None
