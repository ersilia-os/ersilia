import urllib.parse
import requests
import json

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


class CompoundIdentifier(object):
    UNPROCESSABLE_INPUT = "UNPROCESSABLE_INPUT"
    
    def __init__(self, local=True):
        if local:
            self.Chem = Chem
        else:
            self.Chem = None
        self.unichem = unichem
        self.default_type = "smiles"
        self.input_header_synonyms = set(["smiles", "input"])
        self.key_header_synonyms = set(["inchikey", "key"])

    def is_input_header(self, h):
        if h.lower() in self.input_header_synonyms:
            return True
        else:
            return False

    def is_key_header(self, h):
        if h.lower() in self.key_header_synonyms:
            return True
        else:
            return False

    def _is_smiles(self, text):
        if self.Chem is None:
            if self._pubchem_smiles_to_inchikey(text) is not None:
                return True
            else:
                return False
        else:
            mol = self.Chem.MolFromSmiles(text)
            if mol is None:
                return False
            else:
                return True

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
        if text is None:
            return self.default_type
        if self._is_inchikey(text):
            return "inchikey"
        if self._is_smiles(text):
            return "smiles"
        return "name"

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
    def _nci_smiles_to_inchikey(smiles):
        identifier = urllib.parse.quote(smiles)
        url = "https://cactus.nci.nih.gov/chemical/structure/{0}/stdinchikey".format(
            identifier
        )
        req = requests.get(url)
        if req.status_code != 200:
            return None
        return req.text.split("=")[1]

    @staticmethod
    def _pubchem_smiles_to_inchikey(smiles):
        identifier = urllib.parse.quote(smiles)
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{0}/property/InChIKey/json".format(
            identifier
        )
        req = requests.get(url)
        if req.status_code != 200:
            return None
        return json.loads(req.text)["PropertyTable"]["Properties"][0]["InChIKey"]

    @staticmethod
    def chemical_identifier_resolver(identifier):
        """Returns SMILES string of a given identifier, using NCI tool"""
        identifier = urllib.parse.quote(identifier)
        url = "https://cactus.nci.nih.gov/chemical/structure/{0}/smiles".format(
            identifier
        )
        req = requests.get(url)
        if req.status_code != 200:
            return None
        return req.text

    def encode(self, smiles):
        """Get InChIKey of compound based on SMILES string"""
        if smiles is None or not smiles.strip():
            return self.UNPROCESSABLE_INPUT
        
        if self.Chem is None:
            inchikey = self._pubchem_smiles_to_inchikey(smiles) or self._nci_smiles_to_inchikey(smiles)
        else:
            try:
                mol = self.Chem.MolFromSmiles(smiles)
                if mol is None:
                    return self.UNPROCESSABLE_INPUT
                inchi = self.Chem.rdinchi.MolToInchi(mol)[0]
                inchikey = self.Chem.rdinchi.InchiToInchiKey(inchi)
            except:
                inchikey = self.UNPROCESSABLE_INPUT
        
        return inchikey or self.UNPROCESSABLE_INPUT

    
    def process_input(self, smiles_list):
        """
         UNPROCESSABLE_INPUT as both the input and key.
        """
        results = []
        for smiles in smiles_list:
            inchikey = self.encode(smiles)
            if inchikey == self.UNPROCESSABLE_INPUT:
                results.append({
                    'input': self.UNPROCESSABLE_INPUT,
                    'key': self.UNPROCESSABLE_INPUT
                })
            else:
                results.append({
                    'input': smiles,
                    'key': inchikey
                })
        return results

Identifier = CompoundIdentifier
