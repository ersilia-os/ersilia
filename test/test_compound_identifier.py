import pytest
from ersilia.utils.identifiers.compound import CompoundIdentifier

@pytest.fixture
def compound_identifier():
    """Fixture to initialize a CompoundIdentifier instance."""
    return CompoundIdentifier()

def test_is_input_header_positive(compound_identifier):
    """Test that valid input headers return True."""
    valid_headers = ["smiles", "input", "SMILES", "INPUT"]
    for header in valid_headers:
        assert compound_identifier.is_input_header(header) is True
