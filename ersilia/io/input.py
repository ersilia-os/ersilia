from ..utils.identifiers import MoleculeIdentifier


class ChemicalInput(object):

    def __init__(self):
        self.mi = MoleculeIdentifier()

    def single(self, inp):
        inp_type = self.mi.guess_type(inp)
        if inp_type == "smiles":
            smi = inp
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
