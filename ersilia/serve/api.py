import os
import requests
import json
import tempfile

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
        self.logger.debug(
            "API {0}:{1} initialized at URL {2}".format(model_id, api_name, url)
        )

    def _post(self, input, output):
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
                return [{"output": output}]
        else:
            self.logger.error(response)
            return None

    def post(self, input, output, batch_size):
        self.logger.debug("Posting to {0}".format(self.api_name))
        self.logger.debug("Batch size {0}".format(batch_size))
        if output is not None:
            self.logger.debug("Output is a file")
            tmp_folder = tempfile.mkdtemp()
            fmt = output.split(".")[-1]
            i = 0
            for input in self.input_adapter.adapt(input, batch_size=batch_size):
                subfile = os.path.join(tmp_folder, "{chunk-0}.{1}".format(i, fmt))
                self._post(input, subfile)
                subfiles += [subfile]
                i += 1
            self.output_adapter.merge(subfiles, output)
            for o in [output]:
                yield o
        else:
            self.logger.debug("Output is not a file")
            for input in self.input_adapter.adapt(input, batch_size=batch_size):
                result = json.loads(self._post(input, output))
                for r in result:
                    yield r
