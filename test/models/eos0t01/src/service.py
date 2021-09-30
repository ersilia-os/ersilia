import random

from typing import List

from bentoml import BentoService, api, artifacts
from bentoml.types import JsonSerializable
from bentoml.adapters import JsonInput
from bentoml.service.artifacts.common import JSONArtifact


@artifacts([JSONArtifact("model")])
class Service(BentoService):
    @api(input=JsonInput(), batch=True)
    def invert(self, input: List[JsonSerializable]):
        input = input[0]
        output = []
        for inp in input:
            inp = inp["input"]
            output += [{"inverted": inp[::-1]}]
        output = {"result": output,
                  "meta": None}
        return [output]

    @api(input=JsonInput(), batch=True)
    def shuffle(self, input: List[JsonSerializable]):
        input = input[0]
        output = []
        for inp in input:
            inp = inp["input"]
            linp = list(inp)
            random.shuffle(linp)
            output += [{"shuffled": "".join(linp)}]
        return [output]
