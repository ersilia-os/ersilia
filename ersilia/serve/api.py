import os
import requests
import json
import tempfile

from ..io.input import GenericInputAdapter
from ..io.output import GenericOutputAdapter
from .. import logger
from .schema import ApiSchema


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
        try:
            self._empty_output = ApiSchema(
                model_id=model_id, config_json=None
            ).empty_output_by_api(api_name=api_name)
        except:
            self.logger.info("No empty output available")
            self._empty_output = None

    def __result_returner(self, result, output):
        if output is None:
            return result
        else:
            self.output_adapter.adapt(result, output)
            return [{"output": output}]

    def _do_post(self, input, output):
        url = "{0}/{1}".format(self.url, self.api_name)
        response = requests.post(url, json=input)
        if response.status_code == 200:
            result_ = response.json()
            result = []
            for i, o in zip(input, result_):
                result += [{"input": i, "output": o}]
            result = json.dumps(result, indent=4)
            return self.__result_returner(result, output)
        else:
            self.logger.error("Status Code: {0}".format(response.status_code))
            return None

    def _post(self, input, output):
        result = self._do_post(input, output)
        if result is None and self._batch_size > 1:
            self.logger.warning(
                "Batch prediction didn't seem to work. Doing predictions one by one..."
            )
            result = []
            for inp_one in input:
                r = self._do_post([inp_one], output=None)
                if r is None:
                    r = [{"input": inp_one, "output": self._empty_output}]
                else:
                    r = json.loads(r)
                result += r
            result = json.dumps(result, indent=4)
            result = self.__result_returner(result, output)
        return result

    def post(self, input, output, batch_size):
        self._batch_size = batch_size
        self.logger.debug("Posting to {0}".format(self.api_name))
        self.logger.debug("Batch size {0}".format(batch_size))
        if output is not None:
            self.logger.debug("Output is a file")
            tmp_folder = tempfile.mkdtemp()
            fmt = output.split(".")[-1]
            output_base = ".".join(os.path.basename(output).split(".")[:-1])
            i = 0
            subfiles = []
            for input in self.input_adapter.adapt(input, batch_size=batch_size):
                subfile = os.path.join(
                    tmp_folder, "{0}-chunk-{1}.{2}".format(output_base, i, fmt)
                )
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
