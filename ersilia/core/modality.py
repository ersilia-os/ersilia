import tempfile
import os
import json
from .base import ErsiliaBase
from ..utils.terminal import run_command, run_command_check_output
from ..utils.ports import find_free_port


class ServingModalityCommands(ErsiliaBase):

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id

    def current_bentoml(self):
        port = get_free_port()
        cmd = ["bentoml", "serve", self.model_id, "--port", port]



class ServingModalityChecker(ServingModalityCommands):

    def __init__(self, model_id, config_json=None):
        ServingModalityCommands.__init__(self, model_id=model_id, config_json=config_json)

    def current_bentoml(self):
        cmd = ["ersilia", "serve", self.model_id, "--port", port]
        return run_command_check_output(cmd, quiet=True)





class ServingModality(ServingModalityChecker):

    def __init__(self, model_id, config_json=None):
        ServingModalityChecker.__init__(self, model_id=model_id, config_json=config_json)

    def _is_remote(self):
        pass

    def _is_bento(self):
        pass

    def _is_pip(self):
        pass

    def _is_conda(self):
        pass

    def _is_docker(self):
        pass

    def _get_info(self):
        """Get info available from the Bento"""
        tmp_folder = tempfile.mkdtemp()
        tmp_file = os.path.join(tmp_folder, "info.json")
        cmd = "bentoml info {0}:{1} > {2}".format(self.model_id, tag, tmp_file)
        run_command(cmd, quiet=True)
        with open(tmp_file, "r") as f:
            info = json.load(f)
        return info

    def _get_apis(self):
        """Get APIs available for the model, according to the info Bento"""
        info = self._get_info()
        for item in info["apis"]:
            yield item["name"]
