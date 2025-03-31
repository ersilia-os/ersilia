import asyncio
import hashlib
import re

import nest_asyncio

from ...default import UNPROCESSABLE_INPUT

nest_asyncio.apply()


class CompoundIdentifier(object):
    """
    A class to handle compound identification and conversion between different chemical identifiers.

    Parameters
    ----------
    concurrency_limit : int, optional
        The maximum number of concurrent API requests. Default is 10.

    Examples
    --------
    .. code-block:: python

        identifier = CompoundIdentifier()
        smiles = "CCO"
        key = identifier.encode(smiles)

    """

    def __init__(self, concurrency_limit=10):
        self.default_type = "smiles"
        self.input_header_synonyms = set(["smiles", "input"])
        self.key_header_synonyms = set(["key"])
        # The APIs have the worst rate limitation so its better not to increase more than 10
        self.concurrency_limit = concurrency_limit

    def is_input_header(self, h):
        """
        Check if the given header is an input header.

        Parameters
        ----------
        h : str
            The header to check.

        Returns
        -------
        bool
            True if the header is an input header, False otherwise.
        """
        return h.lower() in self.input_header_synonyms

    def is_key_header(self, h):
        """
        Check if the given header is a key header.

        Parameters
        ----------
        h : str
            The header to check.

        Returns
        -------
        bool
            True if the header is a key header, False otherwise.
        """
        return h.lower() in self.key_header_synonyms

    def _is_smiles(self, text):
        if not isinstance(text, str) or not text.strip():
            return False
        # Rough SMILES pattern, allows space for CXSMILES
        SMILES_REGEX = re.compile(
            r"^[A-Za-z0-9@+\-\[\]\(\)=#$:.\\/%,*]+(?:\s[A-Za-z0-9@+\-\[\]\(\)=#$:.\\/%,*]+)*$"
        )
        return bool(SMILES_REGEX.fullmatch(text.strip()))

    def _is_key(self, text):
        if not isinstance(text, str) or not text.strip():
            return False
        KEY_REGEX = re.compile(r"^[A-Za-z0-9]{32}$")
        return bool(KEY_REGEX.fullmatch(text.strip()))

    def convert_smiles_to_checksum(self, smiles):
        """
        Convert a SMILES string to a checksum.

        Parameters
        ----------
        smiles : str
            The SMILES string to convert.

        Returns
        -------
        str
            The converted checksum.
        """
        if smiles is None:
            return None
        return hashlib.md5(smiles.encode()).hexdigest()

    def guess_type(self, text):
        """
        Guess the type of the given text (either 'smiles').

        Parameters
        ----------
        text : str
            The text to guess the type of.

        Returns
        -------
        str
            The guessed type ('smiles', 'UNPROCESSABLE_INPUT').
        """
        if not isinstance(text, str) or not text.strip() or text == UNPROCESSABLE_INPUT:
            return UNPROCESSABLE_INPUT
        if self._is_smiles(text):
            return "smiles"
        return UNPROCESSABLE_INPUT

    async def process_smiles(self, smiles, semaphore, result_list):
        """
        Process a SMILES string asynchronously.

        Parameters
        ----------
        smiles : str
            The SMILES string to process.
        semaphore : asyncio.Semaphore
            The semaphore to limit concurrency.
        session : aiohttp.ClientSession
            The HTTP session for making requests.
        result_list : list
            The list to store results.
        """
        async with semaphore:  # high performance resource manager
            key = self.convert_smiles_to_checksum(smiles)
            result_list.append({"key": key, "input": smiles, "text": smiles})

    async def encode_batch(self, smiles_list):
        """
        Encode a batch of SMILES strings asynchronously.

        Parameters
        ----------
        smiles_list : list
            The list of SMILES strings to encode.

        Returns
        -------
        list
            The list of encoded results.
        """
        result_list = []
        semaphore = asyncio.Semaphore(self.concurrency_limit)
        tasks = []
        for _, smiles in enumerate(smiles_list):
            tasks.append(self.process_smiles(smiles, semaphore, result_list))

        await asyncio.gather(*tasks)
        return result_list

    def encode(self, smiles):
        """
        Get the key of a compound based on its SMILES string.

        Parameters
        ----------
        smiles : str
            The SMILES string of the compound.

        Returns
        -------
        str
            The key of the compound, or 'UNPROCESSABLE_INPUT' if conversion fails.
        """
        if (
            not isinstance(smiles, str)
            or not smiles.strip()
            or smiles == UNPROCESSABLE_INPUT
        ):
            return UNPROCESSABLE_INPUT

        key = self.convert_smiles_to_checksum(smiles)
        return key if key else UNPROCESSABLE_INPUT

    def validate_smiles(self, smiles):
        """
        Validate a SMILES string.

        Parameters
        ----------
        smiles : str
            The SMILES string to validate.

        Returns
        -------
        bool
            True if the SMILES string is valid, False otherwise.
        """
        return self._is_smiles(smiles)


Identifier = CompoundIdentifier
