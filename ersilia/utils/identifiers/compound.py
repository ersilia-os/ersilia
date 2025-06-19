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
        input = "CCO"
        key = identifier.encode(input)

    """

    def __init__(self, concurrency_limit=10):
        self.default_type = "input"
        self.input_header_synonyms = set(["input", "input"])
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

    def _is_input(self, text):
        if not isinstance(text, str) or not text.strip():
            return False
        # Rough input pattern, allows space for CXinput
        input_REGEX = re.compile(
            r"^[A-Za-z0-9@+\-\[\]\(\)=#$:.\\/%,*]+(?:\s[A-Za-z0-9@+\-\[\]\(\)=#$:.\\/%,*]+)*$"
        )
        return bool(input_REGEX.fullmatch(text.strip()))

    def _is_key(self, text):
        if not isinstance(text, str) or not text.strip():
            return False
        KEY_REGEX = re.compile(r"^[A-Za-z0-9]{32}$")
        return bool(KEY_REGEX.fullmatch(text.strip()))

    def convert_input_to_checksum(self, input):
        """
        Convert a input string to a checksum.

        Parameters
        ----------
        input : str
            The input string to convert.

        Returns
        -------
        str
            The converted checksum.
        """
        if input is None:
            return None
        return hashlib.md5(input.encode()).hexdigest()

    def guess_type(self, text):
        """
        Guess the type of the given text (either 'input').

        Parameters
        ----------
        text : str
            The text to guess the type of.

        Returns
        -------
        str
            The guessed type ('input', 'UNPROCESSABLE_INPUT').
        """
        if not isinstance(text, str) or not text.strip() or text == UNPROCESSABLE_INPUT:
            return UNPROCESSABLE_INPUT
        if self._is_input(text):
            return "input"
        return UNPROCESSABLE_INPUT

    async def process_input(self, input, semaphore, result_list):
        """
        Process a input string asynchronously.

        Parameters
        ----------
        input : str
            The input string to process.
        semaphore : asyncio.Semaphore
            The semaphore to limit concurrency.
        session : aiohttp.ClientSession
            The HTTP session for making requests.
        result_list : list
            The list to store results.
        """
        async with semaphore:  # high performance resource manager
            key = self.convert_input_to_checksum(input)
            result_list.append({"key": key, "input": input, "text": input})

    async def encode_batch(self, input_list):
        """
        Encode a batch of input strings asynchronously.

        Parameters
        ----------
        input_list : list
            The list of input strings to encode.

        Returns
        -------
        list
            The list of encoded results.
        """
        result_list = []
        semaphore = asyncio.Semaphore(self.concurrency_limit)
        tasks = []
        for _, input in enumerate(input_list):
            tasks.append(self.process_input(input, semaphore, result_list))

        await asyncio.gather(*tasks)
        return result_list

    def encode(self, input):
        """
        Get the key of a compound based on its input string.

        Parameters
        ----------
        input : str
            The input string of the compound.

        Returns
        -------
        str
            The key of the compound, or 'UNPROCESSABLE_INPUT' if conversion fails.
        """
        if (
            not isinstance(input, str)
            or not input.strip()
            or input == UNPROCESSABLE_INPUT
        ):
            return UNPROCESSABLE_INPUT

        key = self.convert_input_to_checksum(input)
        return key if key else UNPROCESSABLE_INPUT

    def validate_input(self, input):
        """
        Validate a input string.

        Parameters
        ----------
        input : str
            The input string to validate.

        Returns
        -------
        bool
            True if the input string is valid, False otherwise.
        """
        return self._is_input(input)


Identifier = CompoundIdentifier
