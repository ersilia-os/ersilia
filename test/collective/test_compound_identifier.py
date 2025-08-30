from ersilia.default import UNPROCESSABLE_INPUT
import pytest
from ersilia.utils.identifiers.compound import CompoundIdentifier
from unittest.mock import patch


@pytest.fixture
def compound_identifier():
    return CompoundIdentifier()

class TestCompoundIdentifier:

    @pytest.mark.parametrize("header", ["output", "invalid", "InChI", "inchikey", "random"])
    def test_is_input_header_negative(self, compound_identifier, header):
        """Test that invalid input headers return False."""
        assert not compound_identifier.is_input_header(header)

    @pytest.mark.parametrize(
        "header", ["id", "smiles", "inchi", "input", "some_header", "random", "header", ""]
    )
    def test_is_key_header_negative(self, compound_identifier, header):
        assert not compound_identifier.is_key_header(header)

    @pytest.mark.parametrize(
        "smiles, expected", 
        [
            ("CCO", True),
            ("cccc", True),
            ("invalid_smiles", False),
            ("", False)
        ],
    )
    def test_is_smiles(self, compound_identifier, smiles, expected):
        """Test _is_smiles returns True for valid SMILES strings."""
        compound_identifier.Chem = None
        assert compound_identifier._is_smiles(smiles) == expected

    @pytest.mark.parametrize(
        "input, expected",
        [
            ("C[C@@H](O)[C@H](N)C(O)=O", "smiles"),
            ("CCO", "smiles"),
            (None, UNPROCESSABLE_INPUT),
            ("", UNPROCESSABLE_INPUT),
            (" ", UNPROCESSABLE_INPUT),
            ("\n", UNPROCESSABLE_INPUT),
            ("\t", UNPROCESSABLE_INPUT),
            (12345, UNPROCESSABLE_INPUT),
            (3.14, UNPROCESSABLE_INPUT),
            ("𠜎𠜱𡿺𠬠", UNPROCESSABLE_INPUT),        ],
    )
    def test_guess_type(self, compound_identifier, input, expected):
        print(input)
        """Ensure guess_type correctly identifies valid InChIKeys or SMILES strings."""
        assert compound_identifier.guess_type(input) == expected
