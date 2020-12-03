"""Generate identifiers for entities that are relevant to Ersilia, such as compounds or proteins."""

import uuid
import datetime
from hashids import Hashids
from Crypto.Hash import MD5
from datetime import datetime
import random
import string
import urllib.parse
import requests


class IdentifierGenerator(object):

    def __init__(self):
        pass


class LongIdentifier(IdentifierGenerator):

    def __init__(self):
        super().__init__()

    @staticmethod
    def encode():
        """Get UUID code (long identifier)"""
        return str(uuid.uuid4())


class ShortIdentifier(IdentifierGenerator):

    def __init__(self):
        super().__init__()
        self.hashids = Hashids(salt="ersilia is open source", alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890")

    def encode(self):
        """Short identifier based on timestamp"""
        return str(self.hashids.encode(int(datetime.today().timestamp())))


class ModelIdentifier(IdentifierGenerator):

    def __init__(self):
        super().__init__()
        self.letters = string.ascii_lowercase
        self.numbers = "0123456789"

    def encode(self):
        result_str = "eos"
        result_str += random.choice(self.numbers[1:])
        result_str += "".join(random.choice(self.letters+self.numbers) for _ in range(3))
        return result_str

    def generate(self, n):
        ids = set()
        while True:
            if len(ids) < n:
                ids.update([self.encode()])
            else:
                break
        ids = list(ids)
        random.shuffle(ids)
        return ids


class FileIdentifier(IdentifierGenerator):

    def __init__(self, chunk_size=10000):
        super().__init__()
        self.chunk_size = chunk_size

    def encode(self, filename, n=None):
        h = MD5.new()
        with open(filename, "rb") as f:
            while True:
                chunk = f.read(self.chunk_size)
                if len(chunk):
                    h.update(chunk)
                else:
                    break
        return h.hexdigest()[:n]


class MoleculeIdentifier(IdentifierGenerator):

    def __init__(self):
        super().__init__()
        try:
            from rdkit import Chem
            self.Chem = Chem
        except ModuleNotFoundError as err:
            # TODO Logging
            self.Chem = None
        from chembl_webresource_client.unichem import unichem_client as unichem
        self.unichem = unichem

    def _is_smiles(self, text):
        if self.Chem is None:
            # TODO Logging
            pass
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
            return "smiles"
        return "name"

    def unichem_resolver(self, inchikey):
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
        identifier = urllib.parse.quote(identifier)
        url = "https://cactus.nci.nih.gov/chemical/structure/{0}/smiles".format(identifier)
        req = requests.get(url)
        if req.status_code != 200:
            return None
        return req.text

    def encode(self, smiles):
        """Get InChIKey of compound based on SMILES string"""
        if self.Chem is None:
            return None
        mol = self.Chem.MolFromSmiles(smiles)
        if mol is None:
            raise Exception("The SMILES string: %s is not valid or could not be converted to an InChIKey" % smiles)
        inchi = self.Chem.rdinchi.MolToInchi(mol)[0]
        if inchi is None:
            raise Exception("Could not obtain InChI")
        inchikey = self.Chem.rdinchi.InchiToInchiKey(inchi)
        return inchikey


class ProteinIdentifier(IdentifierGenerator):

    def __init__(self):
        super().__init__()
        from Bio.SeqUtils.CheckSum import seguid
        from bioservices.uniprot import UniProt
        self.seguid = seguid
        self.uniprot = UniProt(verbose=False)

    def sequence_from_uniprot(self, uniprot_ac):
        """Returns protein sequence from uniprot identifier"""
        try:
            return self.uniprot.get_fasta_sequence(uniprot_ac)
        except ValueError:
            return None

    @staticmethod
    def protein_identifier_resolver():
        """Returns protein sequence of a given identifier, using """
        pass # TODO

    @staticmethod
    def encode(sequence):
        """Protein seguid checksum based on amino-acid sequence"""
        return str(self.seguid(sequence))
