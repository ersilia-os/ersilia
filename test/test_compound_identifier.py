import pytest
from ersilia.utils.identifiers.compound import CompoundIdentifier

@pytest.fixture
def compound_identifier():
    """Fixture to initialize a CompoundIdentifier instance."""
    return CompoundIdentifier()

def test_is_input_header_positive(compound_identifier):
    """Test that valid input headers return True."""
    assert compound_identifier.is_input_header("smiles") is True
    assert compound_identifier.is_input_header("input") is True
    assert compound_identifier.is_input_header("SMILES") is True
    assert compound_identifier.is_input_header("INPUT") is True
