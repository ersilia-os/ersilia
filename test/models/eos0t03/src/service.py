import random
import json
import collections
from rdkit import Chem
from rdkit.Chem import Descriptors

from bentoml import BentoService, api, artifacts
from bentoml.adapters import JsonInput
from bentoml.service.artifacts.common import JSONArtifact
from bentoml.types import JsonSerializable


@artifacts([JSONArtifact("model")])
class Service(BentoService):
    """The Service class uses BentoService to build multiple inference APIs"""

    @api(input=JsonInput())
    def predict(self, input: JsonSerializable):
        """
        This is a dummy test model.
        It returns the molecular weight of a molecule.
        """
        mol = Chem.MolFromSmiles(input)
        mw = Descriptors.MolWt(mol)
        return json.dumps(mw)
