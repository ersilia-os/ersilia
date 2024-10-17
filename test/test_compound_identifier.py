import pytest
from ersilia.utils.identifiers.compound import CompoundIdentifier

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

@pytest.mark.parametrize("header", ["inchikey", "key", "INCHIKEY", "KEY"])
def test_is_key_header_positive(compound_identifier, header):
    """Test that valid key headers return True."""
    assert compound_identifier.is_key_header(header) is True

@pytest.mark.parametrize("header", ["smiles", "input", "invalid", "123","someHeader"])
def test_is_key_header_negative(compound_identifier, header):
    """Test that invalid key headers return False."""
    assert compound_identifier.is_key_header(header) is False

@pytest.mark.parametrize("header", ["output", "invalid", "InChI", "inchikey"])
def test_is_input_header_negative(compound_identifier, header):
    """Test that invalid input headers return False."""
    assert compound_identifier.is_input_header(header) is False

@pytest.mark.parametrize("header", ["inchikey", "key", "INCHIKEY", "KEY"])
def test_is_key_header_positive(compound_identifier, header):
    """Test that valid key headers return True."""
    assert compound_identifier.is_key_header(header) is True

@pytest.mark.parametrize("header", ["smiles", "input", "invalid", "123","someHeader"])
def test_is_key_header_negative(compound_identifier, header):
    """Test that invalid key headers return False."""
    assert compound_identifier.is_key_header(header) is False