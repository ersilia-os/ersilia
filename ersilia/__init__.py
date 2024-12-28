import os
import warnings

from ._version import __version__
from .core.base import ErsiliaBase
from .core.model import ErsiliaModel
from .core.modelbase import ModelBase
from .default import CONFIG_JSON, EOS, INSTALL_STATUS_FILE, bashrc_cli_snippet
from .utils.config import Config
from .utils.exceptions_utils.throw_ersilia_exception import throw_ersilia_exception
from .utils.logging import logger

# Filter out some warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Disable GPU
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

# Default variables

# Config
if not os.path.exists(os.path.join(EOS, CONFIG_JSON)):
    from .utils.config import Checker

    Checker().config()

# Environmental variables
os.environ["EOS_HOME"] = EOS

# User profile
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

__all__ = [
    "__version__",
    "logger",
    "throw_ersilia_exception",
    "ErsiliaBase",
    "ErsiliaModel",
    "ModelBase",
    "Config",
]
