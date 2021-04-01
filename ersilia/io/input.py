import random
import os
from ..utils.identifiers.molecule import MoleculeIdentifier
from ..utils.identifiers.protein import ProteinIdentifier

MOLECULE_EXAMPLES = "molecules.tsv"
PROTEIN_EXAMPLES = "proteins.tsv"

PATH = os.path.dirname(os.path.realpath(__file__))


class MoleculeInput(object):

    def __init__(self):
        self.mi = MoleculeIdentifier()
        self.example_file = os.path.join(PATH, "examples", MOLECULE_EXAMPLES)

    def random(self):
        with open(self.example_file, "r") as f:
            line = random.choice(f.readlines()).rstrip("\n").split("\t")
        mol = {
            "inchikey": line[0],
            "smiles": line[1],
            "name": line[2]
        }
        return mol

    def single(self, inp):
        inp_type = self.mi.guess_type(inp)
        if inp_type == "smiles":
            smi = inp
        elif inp_type == "inchikey":
            smi = self.mi.unichem_resolver(inp)
        else:
            smi = self.mi.chemical_identifier_resolver(inp)
        return smi

    def multiple(self, filename): # TODO
        R = []
        with open(filename, "r") as f:
            for r in R:
                if not self.mi._is_smiles(r):
                    continue
                else:
                    pass


class ProteinInput(object):

    def __init__(self):
        self.pi = ProteinIdentifier()
        self.example_file = os.path.join(PATH, "examples", PROTEIN_EXAMPLES)

    def random(self):
        with open(self.example_file, "r") as f:
            line = random.choice(f.readlines()).rstrip("\n").split("\t")
        prot = {
            "uniprot_ac": line[0],
            "seq": line[1],
            "name": line[2]
        }
        return prot

    def single(self, inp): # TODO
        pass

    def multiple(self, filename): # TODO
        pass
