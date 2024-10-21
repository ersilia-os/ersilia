import pytest
from ersilia.utils.identifiers.compound import CompoundIdentifier

# testing for _is_smiles positive test when Chem is None.
def test_is_smiles_positive():
    parser = CompoundIdentifier()
    parser.Chem = None
    # Using Dimethyl ether SMILES representation as a valid string
    smiles_text = "COC"

    assert parser._is_smiles(smiles_text) is True