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
    def __init__(self, local=True):
        if local:
            self.Chem = Chem
        else:
            self.Chem = None
        self.unichem = unichem
        self.default_type = "smiles"

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
        if self.Chem is None:
            inchikey = self._pubchem_smiles_to_inchikey(smiles)
            if inchikey is None:
                inchikey = self._nci_smiles_to_inchikey(smiles)
        else:
            try:
                mol = self.Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise Exception(
                        "The SMILES string: %s is not valid or could not be converted to an InChIKey"
                        % smiles
                    )
                inchi = self.Chem.rdinchi.MolToInchi(mol)[0]
                if inchi is None:
                    raise Exception("Could not obtain InChI")
                inchikey = self.Chem.rdinchi.InchiToInchiKey(inchi)
            except:
                inchikey = None
        return inchikey
