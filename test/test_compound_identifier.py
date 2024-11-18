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
    
    
@pytest.mark.parametrize("inchikey", [
    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    "BQJCRHHNABKAKU-KBQPJGBKSA-N",
    "ZJPODVODJYKFHM-UHFFFAOYSA-N"
])
def test_is_inchikey_positive(compound_identifier, inchikey):
    """Test that valid InChIKeys return True."""
    assert compound_identifier._is_inchikey(inchikey) is True


@pytest.mark.parametrize("inchikey", [
    "BSYNRYMUTXBXSQUHFFFAOYSA",        
    "BSYNRYMUTXBXSQ-UHFFFAOYSA-XY", 
    "12345678901234-1234567890-X",
    "BSYNRYMUTXBXSQ_UHFFFAOYSA-N",
    "BSYNRYMUTXBXSQ-UHFFFAOYSA"
])
def test_is_inchikey_negative(compound_identifier, inchikey):
    """Test that invalid InChIKeys return False."""
    assert not compound_identifier._is_inchikey(inchikey)

    
def test_guess_type_with_inchikey(compound_identifier):
    inchikey = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    result = compound_identifier.guess_type(inchikey)
    assert result == "inchikey"
    
    
@patch('ersilia.utils.identifiers.compound.CompoundIdentifier._pubchem_smiles_to_inchikey')
def test_is_smiles_positive_chem_none(mock_pubchem, compound_identifier):
    compound_identifier.Chem = None
    mock_pubchem.return_value = "InChIKey"
    
# Test with a valid SMILES input
    smiles_string = 'CCO' #Ethanol SMILES
    assert compound_identifier._is_smiles(smiles_string) is True
   
    
@patch('requests.get')
async def test_nci_smiles_to_inchikey_positive(mock_get, compound_identifier):
    """Test _nci_smiles_to_inchikey with a mocked positive response."""
    mock_response = mock_get.return_value
    mock_response.status_code = 200
    mock_response.text = "InChIKey=BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    inchikey = await compound_identifier._nci_smiles_to_inchikey(session=None, smiles="CCO")
    assert inchikey == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

  
@patch('requests.get')
async def test_nci_smiles_to_inchikey_negative(mock_get, compound_identifier):
    """Test _nci_smiles_to_inchikey with a mocked negative response."""
    mock_response = mock_get.return_value
    mock_response.status_code = 404  

    inchikey = await compound_identifier._nci_smiles_to_inchikey(session=None, smiles="invalid_smiles")
    assert inchikey is None
