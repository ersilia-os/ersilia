import pytest
from ersilia.utils.identifiers.compound import CompoundIdentifier

@pytest.fixture
def compound_identifier():
    return CompoundIdentifier()

def test_is_key_header_positive(compound_identifier):
    h = compound_identifier
    assert h.is_key_header("key") is True
    assert h.is_key_header("KEY") is True
    assert h.is_key_header("kEY") is True
    assert h.is_key_header("inchikey") is True
    assert h.is_key_header("INCHIKEY") is True
    assert h.is_key_header("inchiKEy") is True
