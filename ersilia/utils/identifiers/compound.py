import asyncio
import nest_asyncio
import aiohttp
import urllib.parse
import requests
from functools import lru_cache
from ..logging import logger

try:
    from chembl_webresource_client.unichem import unichem_client as unichem
except:
    unichem = None
try:
    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
except:
    Chem = None

from ...default import UNPROCESSABLE_INPUT

nest_asyncio.apply()

class CompoundIdentifier(object):
    def __init__(self, local=True, concurrency_limit=10, cache_maxsize=128):
        if local:
            self.Chem = Chem
        else:
            self.Chem = None
        self.unichem = unichem
        self.default_type = "smiles"
        self.input_header_synonyms = set(["smiles", "input"])
        self.key_header_synonyms = set(["inchikey", "key"])
        # The APIs have the worst rate limitation so its better not to increase more than 10
        self.concurrency_limit = concurrency_limit
        self.cache_maxsize = cache_maxsize
        # defining the cache to be dynamic
        self._pubchem_smiles_to_inchikey = lru_cache(maxsize=self.cache_maxsize)(
            self._pubchem_smiles_to_inchikey
        )
        self._nci_smiles_to_inchikey = lru_cache(maxsize=self.cache_maxsize)(
            self._nci_smiles_to_inchikey
        )
        self.chemical_identifier_resolver = lru_cache(maxsize=self.cache_maxsize)(
            self.chemical_identifier_resolver
        )
        self.convert_smiles_to_inchikey_with_rdkit = lru_cache(
            maxsize=self.cache_maxsize
        )(self.convert_smiles_to_inchikey_with_rdkit)

    def is_input_header(self, h):
        return h.lower() in self.input_header_synonyms

    def is_key_header(self, h):
        return h.lower() in self.key_header_synonyms

    def _is_smiles(self, text):
        if not isinstance(text, str) or not text.strip():
            return False  
        if self.Chem is None:
            return asyncio.run(self._process_pubchem_inchikey(text)) is not None
        else:
            mol = self.Chem.MolFromSmiles(text)
            return mol is not None

    async def _process_pubchem_inchikey(self, text):
        async with aiohttp.ClientSession() as session:
            return await self._pubchem_smiles_to_inchikey(session, text)
            
    @staticmethod
    def _is_inchikey(text):
        if len(text) != 27:
            return False
        comp = text.split("-")
        if len(comp) != 3:
            return False
        if len(comp[0]) != 14 or len(comp[1]) != 10 or len(comp[2]) != 1:
            return False
        for c in comp:
            for x in c:
                if not x.isalpha():
                    return False
        return True

    def guess_type(self, text):
        if not isinstance(text, str) or not text.strip() or text == UNPROCESSABLE_INPUT:
            return UNPROCESSABLE_INPUT
        if self._is_inchikey(text):
            return "inchikey"
        if self._is_smiles(text):
            return "smiles"
        return UNPROCESSABLE_INPUT

    def unichem_resolver(self, inchikey):
        if Chem is None or unichem is None:
            return None
        try:
            ret = self.unichem.inchiFromKey(inchikey)
        except:
            return None
        inchi = ret[0]["standardinchi"]
        mol = self.Chem.inchi.MolFromInchi(inchi)
        return self.Chem.MolToSmiles(mol)

    @staticmethod
    def chemical_identifier_resolver(identifier):
        """Returns SMILES string of a given identifier, using NCI tool"""
        if not identifier or not isinstance(identifier, str):
            return UNPROCESSABLE_INPUT
        identifier = urllib.parse.quote(identifier)
        url = "https://cactus.nci.nih.gov/chemical/structure/{0}/smiles".format(
            identifier
        )
        req = requests.get(url)
        if req.status_code != 200:
            return None
        return req.text

    @staticmethod
    async def _pubchem_smiles_to_inchikey(session, smiles):
        """
        Fetch InChIKey for a single SMILES using PubChem API asynchronously.
        The cache is used to store the results of the requests for unique SMILES.
        """
        identifier = urllib.parse.quote(smiles)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{identifier}/property/InChIKey/json"
        try:
            async with session.get(url) as response:
                if response.status != 200:
                    return None
                data = await response.json()
                return data["PropertyTable"]["Properties"][0]["InChIKey"]
        except Exception as e:
            return None

        
    @staticmethod
    async def _nci_smiles_to_inchikey(session, smiles):
        """
        Fetch InChIKey for a single SMILES using NCI asynchronously.
        The cache is used to store the results of the requests for unique SMILES.
        """
        identifier = urllib.parse.quote(smiles)
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{0}/property/InChIKey/json".format(
            identifier
        )
        try:
            async with session.get(url) as response:
                if response.status != 200:
                    return None
                text = await response.text()
                return text.split("=")[1]
        except Exception as e:
            logger.info(f"Failed to fetch InChIKey from NCI for {smiles}: {e}")
            return None
    
    def convert_smiles_to_inchikey_with_rdkit(self, smiles):
        """
        Converts a SMILES string to an InChIKey using RDKit.
        The results are cached to improve performance.
        """
        if not self.Chem:
            return None

        try:
            mol = self.Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES string: {smiles}")
            inchi = self.Chem.rdinchi.MolToInchi(mol)[0]
            inchikey = self.Chem.rdinchi.InchiToInchiKey(inchi)
            return inchikey
        except Exception as e:
            logger.info(f"RDKit failed to convert SMILES {smiles}: {e}")
            return None

    async def process_smiles(self, smiles, semaphore, session, result_list):
        async with semaphore:  # high performance resource manager
            inchikey = self.convert_smiles_to_inchikey_with_rdkit(smiles)

            if inchikey is None:
                inchikey = await self._pubchem_smiles_to_inchikey(session, smiles)
                if inchikey:
                    logger.info("Inchikey converted using PUBCHEM")

            if inchikey is None:
                inchikey = self._nci_smiles_to_inchikey(session, smiles)
                if inchikey:
                    logger.info("Inchikey converted using NCI")

            if inchikey:
                result_list.append({"key": inchikey, "input": smiles, "text": smiles})
            else:
                logger.info(f"No InChIKey found for SMILES {smiles}. Skipping.")

    async def encode_batch(self, smiles_list):
        result_list = []
        semaphore = asyncio.Semaphore(self.concurrency_limit)
        async with aiohttp.ClientSession() as session:
            tasks = []
            for _, smiles in enumerate(smiles_list):
                tasks.append(
                    self.process_smiles(smiles, semaphore, session, result_list)
                )

            await asyncio.gather(*tasks)

        return result_list

    def encode(self, smiles):
            """Get InChIKey of compound based on SMILES string"""
            if not isinstance(smiles, str) or not smiles.strip() or smiles == UNPROCESSABLE_INPUT:
                return UNPROCESSABLE_INPUT

            if self.Chem is None:
                async def fetch_inchikeys():
                    async with aiohttp.ClientSession() as session:
                        inchikey = await self._pubchem_smiles_to_inchikey(session, smiles)
                        if inchikey:
                            return inchikey
                        inchikey = await self._nci_smiles_to_inchikey(session, smiles)
                        return inchikey

                inchikey = asyncio.run(fetch_inchikeys())
            else:
                try:
                    mol = self.Chem.MolFromSmiles(smiles)
                    if mol is None:
                        return UNPROCESSABLE_INPUT
                    inchi = self.Chem.rdinchi.MolToInchi(mol)[0]
                    if inchi is None:
                        return UNPROCESSABLE_INPUT
                    inchikey = self.Chem.rdinchi.InchiToInchiKey(inchi)
                except:
                    inchikey = None
            return inchikey if inchikey else UNPROCESSABLE_INPUT
    
    def validate_smiles(self, smiles):
            return smiles.strip() != "" and Chem.MolFromSmiles(smiles) is not None

Identifier = CompoundIdentifier
