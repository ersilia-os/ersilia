from bentoml import BentoService, api, artifacts
from bentoml.adapters import JsonInput
from bentoml.service.artifacts.common import JSONArtifact
import random
import json
import collections
from pysmiles import read_smiles
from bentoml.types import JsonSerializable


@artifacts([JSONArtifact('model')])
class Service(BentoService):
    """The Service class uses BentoService to build multiple inference APIs"""
    @api(input=JsonInput())
    def predict(self, input: JsonSerializable):
        """
        This is a dummy test model.
        It counts atoms in a SMILES string.
        """
        mol = read_smiles(input, explicit_hydrogen=True)
        counts = collections.defaultdict(int)
        for _, atom in mol.nodes(data="element"):
            counts[atom] += 1
        return json.dumps(counts)
