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
        output = {"result": output, "meta": None}
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

    @api(input=JsonInput(), batch=True)
    def tokenize(self, input: List[JsonSerializable]):
        input = input[0]
        output = []
        for inp in input:
            inp = inp["input"]
            inp = inp[:30].ljust(30, " ")
            linp = list(inp)
            output += [{"tokenized": linp}]
        return [output]

    @api(input=JsonInput(), batch=True)
    def first_characters(self, input: List[JsonSerializable]):
        input = input[0]
        output = []
        for inp in input:
            inp = inp["input"]
            inp = inp[:3].ljust(3, " ")
            linp = list(inp)
            output += [{"fchars": linp}]
        output = {"result": output, "meta": {"fchars": ["char1", "char2", "char3"]}}
        return [output]
