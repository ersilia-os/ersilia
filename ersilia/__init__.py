from ._version import get_versions
from .utils.config import Config
from .core.base import ErsiliaBase
from .core.model import ErsiliaModel


__version__ = get_versions()['version']
del get_versions

__all__ = [
    "__version__"
]

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
