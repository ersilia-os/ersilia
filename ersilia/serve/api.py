import requests
import json

from ..io.input import GenericAdapter


class Api(object):
    def __init__(self, model_id, url, api_name, config_json=None):
        self.model_id = model_id
        self.adapter = GenericAdapter(self.model_id, config_json=config_json)
        self.url = url
        self.api_name = api_name

    def post(self, input):
        input = self.adapter.adapt(input)
        response = requests.post("{0}/{1}".format(self.url, self.api_name), json=input)
        if response.status_code == 200:
            output = response.json()
            result = []
            for i, o in zip(input, output):
                result += [{"input": i, "output": o}]
            result = json.dumps(result, indent=4)
            return result
        else:
            return None
