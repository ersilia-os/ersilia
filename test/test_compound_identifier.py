import pytest
from ersilia.utils.identifiers.compound import CompoundIdentifier

@pytest.fixture
def compound_identifier():
    return CompoundIdentifier()

def test_is_input_header_negative(compound_identifier):
    """Test that invalid input headers return False."""
    invalid_headers = ["name", "ID", "structure", "formula", "xyz", "Smile", "inputt", "SMILE", "inputs"]
    for header in invalid_headers:
        assert compound_identifier.is_input_header(header) is False
