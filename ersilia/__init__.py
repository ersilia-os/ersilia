# External imports
# ruff: noqa
import os
from ._version import __version__
import warnings

# Filter out some warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Disable GPU
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

# Default variables
from .default import EOS, CONFIG_JSON, INSTALL_STATUS_FILE

# Logger
from .utils.logging import logger

# Config
if not os.path.exists(os.path.join(EOS, CONFIG_JSON)):
    from .utils.config import Checker

    Checker().config()

# Exceptions
from .utils.exceptions_utils.throw_ersilia_exception import throw_ersilia_exception

# Environmental variables
os.environ["EOS_HOME"] = EOS

# Global imports
from .utils.config import Config
from .core.base import ErsiliaBase
from .core.modelbase import ModelBase
from .core.model import ErsiliaModel

# User profile
from .default import bashrc_cli_snippet

bashrc_cli_snippet(overwrite=False)


# Check status of installs
def check_install_status():
    fn = os.path.join(EOS, INSTALL_STATUS_FILE)
    if not os.path.exists(fn):
        status = None
    else:
        with open(fn, "r") as f:
            status = f.read().strip()
    results = {"install_status_file": fn, "status": status}
    return results


INSTALL_STATUS = check_install_status()["status"]
