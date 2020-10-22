from ._version import get_versions
from .utils.config import Config
from .core.base import ErsiliaBase
from .core.model import ErsiliaModel
from .utils.installers import check_dependencies


__version__ = get_versions()['version']
del get_versions

check_dependencies()

__all__ = [
    "__version__"
]
