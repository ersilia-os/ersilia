
import json
import os

import requests
import yaml

from .. import ErsiliaBase
from ..default import PACK_METHOD_BENTOML, PACK_METHOD_FASTAPI, PACKMETHOD_FILE
from .paths import Paths


class PackMethodResolver(ErsiliaBase):
    """
    This class resolves the pack method for a given model after the fetch resolution.
    Note that this class is different from the TemplateResolver, which resolves the pack method
    before the fetch procedure.
    """
    def __init__(self, model_id=None, config_json=None):
        self.model_id = model_id
        ErsiliaBase.__init__(self, config_json=config_json)

    def resolve_pack_method_tracked_at_fetch(self):
        """
        Read the pack method from available from the fetch procedure (i.e. as resolved by the TemplateResolver)

        Returns
        -------
        The packaging method (fastapi or bentoml)
        """
        path = self._get_model_path(model_id=self.model_id)
        pack_method_file = os.path.join(path, PACKMETHOD_FILE)
        if not os.path.exists(pack_method_file):
            return None


        pass

    def resolve_pack_method_from_github_metadata(self):
        """
        Resolve the packaging method based on metadata available from GitHub

        Returns
        -------
        The packaging method (fastapi or bentoml)
        """
        model_id = self.model_id
        data = None
        root_github_url = "https://raw.githubusercontent.com/ersilia-os/{0}/refs/heads/main/".format(model_id)
        extensions = ["json", "yml"]
        for ext in extensions:
            url = f"{root_github_url}/metadata.{ext}"
            try:
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    if ext == "json":
                        data = json.loads(response.text)
                    elif ext == "yml":
                        data = yaml.safe_load(response.text)
            except requests.RequestException:
                pass
        if data is None:
            return PACK_METHOD_FASTAPI
        else:
            if "Docker Pack Method" in data.keys():
                return data["Docker Pack Method"].lower()
            else:
                return None

    def resolve_pack_method_source(self):
        """
        Resolve the packaging method for a model based on its source files.

        Parameters
        ----------
        model_path : str
            The path to the model directory.

        Returns
        -------
        str or None
            The packaging method if found, otherwise None.
        """
        model_path = self.model_path
        if os.path.exists(os.path.join(model_path, "installs", "install.sh")):
            return PACK_METHOD_FASTAPI
        elif os.path.exists(os.path.join(model_path, "bentoml.yml")):
            return PACK_METHOD_BENTOML
        self.logger.warning("Could not resolve pack method")
        return None

    def resolve_pack_method(self):
        """
        Resolve the packaging method for a model.

        Parameters
        ----------
        model_path : str
            The path to the model directory.

        Returns
        -------
        str
            The packaging method.
        """
        model_path = self.model_path
        service_class_file = os.path.join(model_path, "service_class.txt")
        if not os.path.exists(service_class_file):
            service_class = "pulled_docker"
        else:
            with open(os.path.join(model_path, "service_class.txt"), "r") as f:
                service_class = f.read().strip()
        if service_class == "pulled_docker":
            model_id = Paths().model_id_from_path(model_path)
            return self.resolve_pack_method_from_github_metadata(model_id)
        else:
            return self.resolve_pack_method_source(model_path)
