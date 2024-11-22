
from ersilia.default import UNPROCESSABLE_INPUT
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

@pytest.fixture(params=[True, False], ids=["Chem_None", "Chem_Not_None"])
def compound_identifier(request):
    """Fixture that initializes CompoundIdentifier with or without RDKit."""
    return CompoundIdentifier(local=request.param)

@pytest.mark.parametrize("smiles, expected", [
    ("C", True), 
    ("CCO", True) 
])
def test_is_smiles_positive(compound_identifier, smiles, expected):
    """Test _is_smiles returns True for valid SMILES strings."""
    if compound_identifier.Chem is None:
        assert compound_identifier._is_smiles(smiles) == expected

@pytest.mark.parametrize("smiles, expected", [
    ("invalid_smiles", False),  
    ("", False)  
])
def test_is_smiles_negative(compound_identifier, smiles, expected):
    """Test _is_smiles returns False for invalid or empty SMILES strings."""
    assert compound_identifier._is_smiles(smiles) == expected

@pytest.mark.parametrize("inchikey, expected", [
    ("BQJCRHHNABKAKU-KBQPJGBKSA-N", True), 
])
def test_is_inchikey_positive(inchikey, expected):
    """Test _is_inchikey returns True for valid InChIKey."""
    assert CompoundIdentifier._is_inchikey(inchikey) == expected

@pytest.mark.parametrize("inchikey, expected", [
    ("invalid_inchikey", False),
    ("BQJCRHHNABKAKU-KBQPJGBKSA", False) 
])
def test_is_inchikey_negative(inchikey, expected):
    """Test _is_inchikey returns False for invalid InChIKeys."""
    assert CompoundIdentifier._is_inchikey(inchikey) == expected

@pytest.mark.parametrize("inchikey, expected", [
    ("BQJCRHHNABKAKU-KBQPJGBKSA-N", "inchikey"),  
    ("ABCDEFGHIJKLMN-OPQRSTUVWX-Y", "inchikey"),  
])
def test_guess_type_inchikey(compound_identifier, inchikey, expected):
    """Ensure guess_type correctly identifies valid InChIKeys."""
    result = compound_identifier.guess_type(inchikey)
    assert result == expected, f"Expected 'inchikey', but got '{result}' for input '{inchikey}'"

@pytest.mark.parametrize("smiles, expected", [
    ("C", "smiles"),  
    ("CCO", "smiles"), 
])
def test_guess_type_smiles(compound_identifier, smiles, expected):
    """Ensure guess_type correctly identifies valid SMILES strings."""
    result = compound_identifier.guess_type(smiles)
    assert result == expected, f"Expected 'smiles', but got '{result}' for input '{smiles}'"

@pytest.mark.parametrize("input_data, expected", [
    (None, UNPROCESSABLE_INPUT),  
    (UNPROCESSABLE_INPUT, UNPROCESSABLE_INPUT),  
])
def test_guess_type_unprocessable(compound_identifier, input_data, expected):
    """Ensure guess_type returns UNPROCESSABLE_INPUT for None or unprocessable inputs."""
    result = compound_identifier.guess_type(input_data)
    assert result == expected, f"Expected '{UNPROCESSABLE_INPUT}', but got '{result}'"

@pytest.mark.parametrize("whitespace_input, expected", [
    ("\n", UNPROCESSABLE_INPUT),  
    ("\t", UNPROCESSABLE_INPUT),  
    (" ", UNPROCESSABLE_INPUT),   
])
def test_guess_type_whitespace(compound_identifier, whitespace_input, expected):
    """Ensure guess_type returns UNPROCESSABLE_INPUT for whitespace-only input."""
    result = compound_identifier.guess_type(whitespace_input)
    assert result == expected, f"Expected '{UNPROCESSABLE_INPUT}' for input '{whitespace_input}'"

@pytest.mark.parametrize("non_char_input, expected", [
    (12345, UNPROCESSABLE_INPUT),  
    (3.14, UNPROCESSABLE_INPUT),   
    ("𠜎𠜱𡿺𠬠", UNPROCESSABLE_INPUT),  
])
def test_guess_type_non_character(compound_identifier, non_char_input, expected):
    """Ensure guess_type returns UNPROCESSABLE_INPUT for non-character input."""
    result = compound_identifier.guess_type(non_char_input)
    assert result == expected, f"Expected '{UNPROCESSABLE_INPUT}' for input '{non_char_input}'"
    
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
