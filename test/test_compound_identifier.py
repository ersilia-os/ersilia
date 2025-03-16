from ersilia.default import UNPROCESSABLE_INPUT
import pytest
from ersilia.utils.identifiers.compound import CompoundIdentifier
from unittest.mock import AsyncMock, patch, MagicMock
import aiohttp
import requests


@pytest.fixture
def compound_identifier():
    return CompoundIdentifier(local=False) 


class TestCompoundIdentifier:

    @pytest.mark.parametrize("header", ["smiles", "input", "SMILES", "INPUT"])
    def test_is_input_header_positive(self, compound_identifier, header):
        """Test that valid input headers return True."""
        assert compound_identifier.is_input_header(header)

    @pytest.mark.parametrize("header", ["output", "invalid", "InChI", "inchikey", "random"])
    def test_is_input_header_negative(self, compound_identifier, header):
        """Test that invalid input headers return False."""
        assert not compound_identifier.is_input_header(header)

    @pytest.mark.parametrize("header", ["key", "inchiKey", "KEY", "INCHIKEY"])
    def test_is_key_header_positive(self, compound_identifier, header):
        """Test that valid key headers return True."""
        assert compound_identifier.is_key_header(header)

    @pytest.mark.parametrize(
        "header", ["id", "smiles", "inchi", "input", "some_header", "random", "header", ""]
    )
    def test_is_key_header_negative(self, compound_identifier, header):
        assert not compound_identifier.is_key_header(header)

    @pytest.mark.parametrize(
        "inchikey, expected",
        [
            ("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", True),
            ("random-text-that-is-not-an-inchikey", False),
        ],
    )
    def test_is_inchikey(self, compound_identifier, inchikey, expected):
        """Test that valid InChIKeys return True."""
        assert compound_identifier._is_inchikey(inchikey) == expected

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
            ("BQJCRHHNABKAKU-KBQPJGBKSA-N", "inchikey"),
            ("ABCDEFGHIJKLMN-OPQRSTUVWX-Y", "inchikey"),
            ("C", "smiles"),
            ("CCO", "smiles"),
            (None, UNPROCESSABLE_INPUT),
            ("", UNPROCESSABLE_INPUT),
            (" ", UNPROCESSABLE_INPUT),
            ("\n", UNPROCESSABLE_INPUT),
            ("\t", UNPROCESSABLE_INPUT),
            (12345, UNPROCESSABLE_INPUT),
            (3.14, UNPROCESSABLE_INPUT),
            ("𠜎𠜱𡿺𠬠", UNPROCESSABLE_INPUT),
            ("random text", UNPROCESSABLE_INPUT)
        ],
    )
    def test_guess_type(self, compound_identifier, input, expected):
        """Ensure guess_type correctly identifies valid InChIKeys or SMILES strings."""
        assert compound_identifier.guess_type(input) == expected

    @pytest.mark.asyncio
    @patch("aiohttp.ClientSession.get")
    async def test_nci_smiles_to_inchikey_positive(self, mock_get, compound_identifier):
        """Test NCI with correct mock"""
        mock_response = MagicMock()
        mock_response.status = 200
        mock_response.text = AsyncMock(return_value="BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        mock_get.return_value.__aenter__.return_value = mock_response

        async with aiohttp.ClientSession() as session:
            inchikey = await compound_identifier._nci_smiles_to_inchikey(session, "CCO")
        assert inchikey == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    @pytest.mark.asyncio
    @patch("requests.get")
    async def test_nci_smiles_to_inchikey_negative(self, mock_get, compound_identifier):
        """Test _nci_smiles_to_inchikey with a mocked negative response."""
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response

        inchikey = await compound_identifier._nci_smiles_to_inchikey(
            session=None, smiles="invalid_smiles"
        )
        assert inchikey is None

    @pytest.mark.asyncio
    @patch("aiohttp.ClientSession.get")
    async def test_pubchem_smiles_to_inchikey_positive(self, mock_get, compound_identifier):
        """Test PubChem with correct mock"""
        mock_response = MagicMock()
        mock_response.status = 200
        mock_response.json = AsyncMock(return_value={
            "PropertyTable": {
                "Properties": [{"InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"}]
            }
        })
        mock_get.return_value.__aenter__.return_value = mock_response

        async with aiohttp.ClientSession() as session:
            inchikey = await compound_identifier._pubchem_smiles_to_inchikey(session, "CCO")
        assert inchikey == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


    @pytest.mark.asyncio
    @patch("requests.get")
    async def test_pubchem_smiles_to_inchikey_negative(self, mock_get, compound_identifier):
        """Test _pubchem_smiles_to_inchikey with a mocked positive response."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"}
        mock_get.return_value = mock_response

        inchikey = await compound_identifier._pubchem_smiles_to_inchikey(
            session=None, smiles="invalid_smiles"
        )
        assert inchikey is None

    @patch("requests.get")
    def test_chemical_identifier_resolver_valid(self, mock_get, compound_identifier):
        """Test chemical_identifier_resolver with valid input"""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "CC(=O)Oc1ccccc1C(O)=O"
        mock_get.return_value = mock_response

        smiles = compound_identifier.chemical_identifier_resolver("aspirin")
        assert smiles == "CC(=O)Oc1ccccc1C(O)=O"

    @patch("requests.get")
    def test_chemical_identifier_resolver_invalid(self, mock_get, compound_identifier):
        """Test chemical_identifier_resolver with incorrect input (mocked 404 response)"""
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response

        result = compound_identifier.chemical_identifier_resolver("someincorrectinput")
        
        assert result is None

    

