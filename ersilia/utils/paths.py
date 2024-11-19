from dataclasses import dataclass, asdict
import re
import os
import json
from typing import List, Optional
import yaml
from pathlib import Path
from ersilia import logger
from .docker import resolve_pack_method_docker
from ..default import PACK_METHOD_BENTOML, PACK_METHOD_FASTAPI, METADATA_JSON_FILE, METADATA_YAML_FILE

MODELS_DEVEL_DIRNAME = "models"


class Paths(object):
    def __init__(self):
        self.essentials = ["setup.py", "README.md", "CODE_OF_CONDUCT.md"]

    @staticmethod
    def _eos_regex():
        return re.compile(r"eos[0-9][a-z0-9]{3}")

    @staticmethod
    def home():
        """Get home directory"""
        return os.path.abspath(str(Path.home()))

    def model_id_from_path(self, path):
        """Guess model identifier based on the path"""
        regex = self._eos_regex()
        path = os.path.abspath(path)
        model_ids = sorted(set(regex.findall(path)))
        if len(model_ids) == 1:
            return model_ids[0]
        else:
            return None

    def org_development_path(self):
        """Guess generic development path"""
        path = self.ersilia_development_path()
        if path is None:
            return
        else:
            return os.path.split(path)[0]

    def ersilia_development_path(self):
        """Try to guess the package development path in the local computer"""
        path = os.path.dirname(__file__)
        for _ in range(2):
            path = os.path.split(path)[0]
        for essential in self.essentials:
            if not os.path.exists(os.path.join(path, essential)):
                return None
        return path

    @staticmethod
    def exists(path):
        if path is None:
            return False
        if os.path.exists(path):
            return True
        else:
            return False

@dataclass(init=True)
class Metadata:
    Identifier: str
    Slug: str
    Title: str
    Description: str
    Mode: str
    Input: List[str]
    InputShape: str
    Task: List[str]
    Output: List[str]
    OutputType: List[str]
    OutputShape: str
    Interpretation: str
    Tag: List[str]
    Publication: str
    SourceCode: str
    License: str
    DockerHub: Optional[str] = None
    DockerArchitecture: Optional[List[str]] = None
    S3: Optional[str] = None
    Status: Optional[str] = None
    Contributor: Optional[str] = None

    def __post_init__(self):
        # Dataclasses do not implicitly perform type checks
        # This is a workaround to ensure that the fields are of the correct type
        for field_name, field_def in self.__dataclass_fields__.items():
            field_value = getattr(self, field_name)
            if field_def.type == List[str] and not isinstance(field_value, list):
                try:
                    setattr(self, field_name, [field_value])
                except TypeError:
                    raise TypeError(f"Field {field_name} has the wrong type")

class ErsiliaMetadataLoader(yaml.SafeLoader):
    """
    YAML loader for Ersilia metadata
    Mainly this is needed so we don't directly modify the SafeLoader class,
    which would break the YAML parsing for other parts of the code
    """
    pass

def metadata_constructor(loader, node):
    _fields = loader.construct_mapping(node)
    fields = {}
    for key in _fields.keys():
        key_split = key.split()
        if len(key_split) > 1:
            fields["".join(key_split)] = _fields[key]
        else:
            fields[key] = _fields[key]
    return Metadata(**fields)

ErsiliaMetadataLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, metadata_constructor
)

def resolve_pack_method_source(model_path):
    if os.path.exists(os.path.join(model_path, "installs", "install.sh")):
        return PACK_METHOD_FASTAPI
    elif os.path.exists(os.path.join(model_path, "bentoml.yml")):
        return PACK_METHOD_BENTOML
    logger.warning("Could not resolve pack method")
    return None

def resolve_pack_method(model_path):
    with open(os.path.join(model_path, "service_class.txt"), "r") as f:
        service_class = f.read().strip()
    if service_class == "pulled_docker":
        model_id = Paths().model_id_from_path(model_path)
        return resolve_pack_method_docker(model_id)
    else:
        return resolve_pack_method_source(model_path)
    
def get_metadata_from_base_dir(path):
    if os.path.exists(os.path.join(path, METADATA_JSON_FILE)):
        with open(os.path.join(path, METADATA_JSON_FILE), "r") as f:
            metadata = json.load(f)
    elif os.path.exists(os.path.join(path, METADATA_YAML_FILE)):
        with open(os.path.join(path, METADATA_YAML_FILE), "r") as f:
            try:
                metadata = asdict(yaml.load(f, Loader=ErsiliaMetadataLoader))
            except TypeError:
                return 
    else:
        raise FileNotFoundError("Metadata file not found")
    return metadata

