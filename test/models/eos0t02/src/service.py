import random
import json
import collections
from pysmiles import read_smiles
from typing import List

from bentoml import BentoService, api, artifacts
from bentoml.adapters import JsonInput
from bentoml.service.artifacts.common import JSONArtifact
from bentoml.types import JsonSerializable


@artifacts([JSONArtifact("model")])
class Service(BentoService):
    @api(input=JsonInput(), batch=True)
    def predict(self, input: List[JsonSerializable]):
        """
        This is a dummy test model.
        It counts atoms in a SMILES string.
        """
        input = input[0]
        output = []
        for inp in input:
            mol = read_smiles(inp["input"], explicit_hydrogen=True)
            counts = collections.defaultdict(int)
            for _, atom in mol.nodes(data="element"):
                counts[atom] += 1
            output += [{"atoms": counts}]
        return [output]
