# Version
from ._version import __version__
del _version

# External imports
import os
import click

# Default variables
from .default import EOS, CONFIG_JSON, INSTALL_STATUS_FILE
if not os.path.exists(os.path.join(EOS, CONFIG_JSON)):
    from .utils.config import Checker
    Checker().config()

# Environmental variables
os.environ["EOS_HOME"] = EOS

# Global imports
from .utils.config import Config
from .core.base import ErsiliaBase
from .core.model import ErsiliaModel

# Clean version
from ._clean_static_version import version
script_path = os.path.dirname(os.path.abspath(__file__))
clean_version_file = os.path.join(script_path, "_clean_static_version.py")
if __version__[:7] == "unknown":
    __version__ = version
else:
    ver = __version__.split(".")
    ver = "{0}.{1}.{2}".format(ver[0], ver[1], ver[2].split("+")[0])
    if ver != version:
        with open(clean_version_file, "w") as f:
            f.write('version = "{0}"'.format(ver))
    __version__ = ver

# Check status of installs
def check_install_status():
    fn = os.path.join(EOS, INSTALL_STATUS_FILE)
    if not os.path.exists(fn):
        status = None
    else:
        with open(fn, "r") as f:
            status = f.read().strip()
    results = {
        "install_status_file": fn,
        "status": status
    }
    return results


INSTALL_STATUS = check_install_status()["status"]


__all__ = [
    "__version__"
]
