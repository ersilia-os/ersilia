import requests
import json

from ..io.input import GenericInputAdapter
from ..io.output import GenericOutputAdapter
from .. import logger


class Api(object):
    def __init__(self, model_id, url, api_name, config_json=None):
        self.model_id = model_id
        self.input_adapter = GenericInputAdapter(self.model_id, config_json=config_json)
        self.output_adapter = GenericOutputAdapter(config_json=config_json)
        if url[-1] == "/":
            self.url = url[:-1]
        else:
            self.url = url
        self.api_name = api_name
        self.logger = logger
        self.logger.debug("API {0}:{1} initialized at URL {2}".format(model_id, api_name, url))

    def post(self, input, output=None):
        logger.debug("Posting")
        input = self.input_adapter.adapt(input)
        url = "{0}/{1}".format(self.url, self.api_name)
        self.logger.debug(url)
        self.logger.debug(input)
        response = requests.post(url, json=input)
        self.logger.debug("{0}".format(response.status_code))
        self.logger.debug(response.text)
        if response.status_code == 200:
            result_ = response.json()
            result = []
            for i, o in zip(input, result_):
                result += [{"input": i, "output": o}]
            result = json.dumps(result, indent=4)
            if output is None:
                return result
            else:
                self.output_adapter.adapt(result, output)
                return {"output": output}
        else:
            self.logger.error(response)
            return None
