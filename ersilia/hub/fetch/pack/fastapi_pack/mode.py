import json
import os

from ..... import ErsiliaBase
from .....default import MODEL_CONFIG_FILENAME
from .....utils.system import SystemChecker

AVAILABLE_MODES = ["system", "conda"]


class PackModeDecision(ErsiliaBase):
    """
    Class to decide the mode of operation for a model based on its configuration and the environment.

    Parameters
    ----------
    model_id : str
        The identifier of the model.
    config_json : dict
        The configuration JSON for the model.

    Examples
    --------
    .. code-block:: python

        pmd = PackModeDecision(
            model_id="model123", config_json={}
        )
        mode = pmd.decide()
    """

    def __init__(self, model_id: str, config_json: dict):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def decide_from_config_file_if_available(self) -> str:
        """
        Decide the mode from the configuration file if available.

        Returns
        -------
        str
            The mode specified in the configuration file, or None if not specified.

        Raises
        ------
        Exception
            If the default_mode in the configuration file is not one of the AVAILABLE_MODES.
        """
        folder = self._model_path(self.model_id)
        if not os.path.exists(os.path.join(folder, MODEL_CONFIG_FILENAME)):
            return None
        with open(os.path.join(folder, MODEL_CONFIG_FILENAME), "r") as f:
            model_config = json.load(f)
        if "default_mode" in model_config:
            default_mode = model_config["default_mode"]
            if default_mode not in AVAILABLE_MODES:
                raise Exception(
                    "The model default_mode specified in the config.json file of the model repo is not correct. It should be one of {0}".format(
                        " ".join(AVAILABLE_MODES)
                    )
                )
            else:
                return default_mode
        return None

    def decide(self) -> str:
        """
        Decide the mode of operation for the model.

        Returns
        -------
        str
            The mode of operation, either "system" or "conda".
        """
        sc = SystemChecker()
        if sc.is_github_action():
            self.logger.debug(
                "Code is being run inside a GitHub Actions workflow. Use conda as a by-default mode."
            )
            return "conda"
        mode = self.decide_from_config_file_if_available()
        if mode is not None:
            self.logger.debug("Mode is already specified in the model repository")
            self.logger.debug("Mode: {0}".format(mode))
            return mode
        self.logger.debug(
            "Check if model can be run with vanilla (system) code. This is a default option when inside a docker container."
        )
        if sc.is_inside_docker():
            return "system"
        else:
            return "conda"
