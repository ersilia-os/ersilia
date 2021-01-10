from ._version import __version__
del _version

# Global imports
from .utils.config import Config
from .core.base import ErsiliaBase
from .core.model import ErsiliaModel
from .utils.installers import check_dependencies

# Clean version
import os
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

# Check dependencies
check_dependencies()

__all__ = [
    "__version__"
]
