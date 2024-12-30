from dataclasses import dataclass, asdict
import re
import os
import json
from typing import List, Optional
import yaml
from pathlib import Path
from ersilia import logger
from .docker import resolve_pack_method_docker
from ..default import (
    PACK_METHOD_BENTOML,
    PACK_METHOD_FASTAPI,
    METADATA_JSON_FILE,
    METADATA_YAML_FILE,
)

MODELS_DEVEL_DIRNAME = "models"


class Paths(object):
    """
    A class to handle various path-related operations in Ersilia.

    Methods
    -------
    model_id_from_path(path)
        Guess model identifier based on the path.
    org_development_path()
        Guess generic development path.
    ersilia_development_path()
        Try to guess the package development path in the local computer.
    exists(path)
        Check if a path exists.
    """
    def __init__(self):
        self.essentials = ["setup.py", "README.md", "CODE_OF_CONDUCT.md"]

    @staticmethod
    def _eos_regex():
        return re.compile(r"eos[0-9][a-z0-9]{3}")

    @staticmethod
    def home():
        """
        Get home directory.

        Returns
        -------
        str
            The home directory path.
        """
        return os.path.abspath(str(Path.home()))

    def model_id_from_path(self, path):
        """
        Guess model identifier based on the path.

        Parameters
        ----------
        path : str
            The path to guess the model identifier from.

        Returns
        -------
        str or None
            The model identifier if found, otherwise None.
        """
        regex = self._eos_regex()
        path = os.path.abspath(path)
        model_ids = sorted(set(regex.findall(path)))
        if len(model_ids) == 1:
            return model_ids[0]
        else:
            return None

    def org_development_path(self):
        """
        Guess generic development path.

        Returns
        -------
        str or None
            The generic development path if found, otherwise None.
        """
        path = self.ersilia_development_path()
        if path is None:
            return
        else:
            return os.path.split(path)[0]

    def ersilia_development_path(self):
        """
        Try to guess the package development path in the local computer.

        Returns
        -------
        str or None
            The package development path if found, otherwise None.
        """
        path = os.path.dirname(__file__)
        for _ in range(2):
            path = os.path.split(path)[0]
        for essential in self.essentials:
            if not os.path.exists(os.path.join(path, essential)):
                return None
        return path

    @staticmethod
    def exists(path):
        """
        Check if a path exists.

        Parameters
        ----------
        path : str
            The path to check.

        Returns
        -------
        bool
            True if the path exists, False otherwise.
        """
        if path is None:
            return False
        if os.path.exists(path):
            return True
        else:
            return False


@dataclass(init=True)
class Metadata:
    """
    A dataclass to represent metadata for Ersilia models.

    Attributes
    ----------
    Identifier : str
        The model identifier.
    Slug : str
        The model slug.
    Title : str
        The model title.
    Description : str
        The model description.
    Mode : str
        The model mode.
    Input : List[str]
        The model input types.
    InputShape : str
        The shape of the model input.
    Task : List[str]
        The tasks the model performs.
    Output : List[str]
        The model output types.
    OutputType : List[str]
        The types of the model output.
    OutputShape : str
        The shape of the model output.
    Interpretation : str
        The interpretation of the model output.
    Tag : List[str]
        The tags associated with the model.
    Publication : str
        The publication associated with the model.
    SourceCode : str
        The source code repository for the model.
    License : str
        The license for the model.
    DockerHub : Optional[str], optional
        The DockerHub repository for the model. Default is None.
    DockerArchitecture : Optional[List[str]], optional
        The Docker architectures supported by the model. Default is None.
    S3 : Optional[str], optional
        The S3 bucket for the model. Default is None.
    Status : Optional[str], optional
        The status of the model. Default is None.
    Contributor : Optional[str], optional
        The contributor of the model. Default is None.
    """
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
    if os.path.exists(os.path.join(model_path, "installs", "install.sh")):
        return PACK_METHOD_FASTAPI
    elif os.path.exists(os.path.join(model_path, "bentoml.yml")):
        return PACK_METHOD_BENTOML
    logger.warning("Could not resolve pack method")
    return None


def resolve_pack_method(model_path):
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
    with open(os.path.join(model_path, "service_class.txt"), "r") as f:
        service_class = f.read().strip()
    if service_class == "pulled_docker":
        model_id = Paths().model_id_from_path(model_path)
        return resolve_pack_method_docker(model_id)
    else:
        return resolve_pack_method_source(model_path)


def get_metadata_from_base_dir(path):
    """
    Get metadata from the base directory of a model.

    Parameters
    ----------
    path : str
        The path to the base directory.

    Returns
    -------
    dict
        The metadata dictionary.

    Raises
    ------
    FileNotFoundError
        If the metadata file is not found.
    """
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
