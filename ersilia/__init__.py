from ._version import __version__
del _version

from .utils.config import Config
from .core.base import ErsiliaBase
from .core.model import ErsiliaModel
from .utils.installers import check_dependencies


check_dependencies()

__all__ = [
    "__version__"
]
