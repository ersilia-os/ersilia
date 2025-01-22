from typing import List

from bentoml import BentoService, api, artifacts
from bentoml.adapters import JsonInput
from bentoml.service.artifacts.common import JSONArtifact
from bentoml.types import JsonSerializable
from rdkit import Chem
from rdkit.Chem import Descriptors


@artifacts([JSONArtifact("model")])
class Service(BentoService):
    @api(input=JsonInput(), batch=True)
    def predict(self, input: List[JsonSerializable]):
        """
        This is a dummy test model.
        It returns the molecular weight of a molecule.
        """
        input = input[0]
        output = []
        for inp in input:
            mol = Chem.MolFromSmiles(inp["input"])
            mw = Descriptors.MolWt(mol)
            output += [{"mw": mw}]
        return [output]
