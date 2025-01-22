import importlib
import os

try:
    from isaura.default import REPOSITORY_PATH as ISAURA_REPOSITORY_PATH
except:
    ISAURA_REPOSITORY_PATH = None

from .. import ErsiliaBase


class LakeBase(ErsiliaBase):
    """
    Base class for managing the lake directory.

    Parameters
    ----------
    config_json : dict
        Configuration settings in JSON format.

    Attributes
    ----------
    lake_dir : str or None
        Absolute path to the lake directory if ISAURA_REPOSITORY_PATH is set, otherwise None.
    """

    def __init__(self, config_json: dict):
        ErsiliaBase.__init__(self, config_json=config_json)
        if ISAURA_REPOSITORY_PATH is not None:
            self.lake_dir = os.path.abspath(ISAURA_REPOSITORY_PATH)
        else:
            self.lake_dir = None

    def is_installed(self) -> bool:
        """
        Check if the 'isaura' package is installed.

        Returns
        -------
        bool
            True if 'isaura' is installed, False otherwise.

        Warns
        -----
        ModuleNotFoundError
            If 'isaura' is not installed, a warning is logged.
        """
        try:
            importlib.util.find_spec("isaura")

            return True
        except ModuleNotFoundError:
            self.logger.warning(
                "Lake manager 'isaura' is not installed! We strongly recommend installing it to store calculations persistently"
            )
            return False
