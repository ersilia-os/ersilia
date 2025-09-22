from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

MW_BINS = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500]
LOGP_BINS = [-1, 0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5]

print(len(MW_BINS), len(LOGP_BINS))


def tranche_coordinates(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)

    for i, edge in enumerate(MW_BINS):
        if mw <= edge:
            col = i
            break
    else:
        col = len(MW_BINS)

    for j, edge in enumerate(LOGP_BINS):
        if logp <= edge:
            row = j
            break
    else:
        row = len(LOGP_BINS)

    return (row, col)


# Example usage:
if __name__ == "__main__":
    smiles = "CC1CCC2C(C1)CCC3C2CC=C4C3(CCC(=O)C5=CC(=C(C=C5C(=O)O4)C)C)C(=O)C=CC6=CC(=C(C=C6C(=O)C))O"
    row, col = tranche_coordinates(smiles)
    print(f"SMILES: {smiles}, Tranche Coordinates: (Row: {row}, Column: {col})")