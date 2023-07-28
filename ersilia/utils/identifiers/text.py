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


class TextIdentifier(object):
    def __init__(self, local=True):
        if local:
            self.Chem = Chem
        else:
            self.Chem = None
        self.unichem = unichem

    def _is_iupac_name(self, text):
        if self._pubchem_iupac_name_to_inchikey(text) is not None:
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
        if self._is_inchikey(text):
            return "inchikey"
        if self._is_smiles(text):
            return "smile"
        if self._is_iupac_name(text):
            return "iupac_name"
        return "name"

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
    def _nci_iupac_name_to_inchikey(iupac_name):
        """Returns inchi key string of a given identifier, using NCI tool"""
        identifier = urllib.parse.quote(iupac_name)
        url = "https://cactus.nci.nih.gov/chemical/structure/{0}/stdinchikey".format(
            identifier
        )
        req = requests.get(url)
        if req.status_code != 200:
            return None
        return req.text.split("=")[1]

    @staticmethod
    def _pubchem_iupac_name_to_inchikey(iupac_name):
        """Returns inchi key string of a given identifier, using PubChem"""
        identifier = urllib.parse.quote(iupac_name)
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{0}/property/InChIKey/json".format(
            identifier
        )
        req = requests.get(url)
        if req.status_code != 200:
            return None
        return json.loads(req.text)["PropertyTable"]["Properties"][0]["InChIKey"]

    @staticmethod
    def chemical_identifier_resolver(identifier):
        """Returns iupac name string of a given identifier, using NCI tool"""
        identifier = urllib.parse.quote(identifier)
        url = "https://cactus.nci.nih.gov/chemical/structure/{0}/iupac_name".format(
            identifier
        )
        req = requests.get(url)
        if req.status_code != 200:
            return None
        return req.text

    def encode(self, iupac_name):
        """Get InChIKey of compound based on iupac name string"""
        inchikey = self._pubchem_iupac_name_to_inchikey(iupac_name)
        if inchikey is None:
            inchikey = self._nci_iupac_name_to_inchikey(iupac_name)
        return inchikey
