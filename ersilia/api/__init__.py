"""
The Ersilia Python API: ``Model`` and ``Catalog`` mirror the CLI commands.

Every error is raised as an ``ErsiliaError`` subclass.
"""

from ..utils.exceptions_utils.exceptions import ErsiliaError
from .create_api import Catalog, Model

__all__ = ["Catalog", "ErsiliaError", "Model"]
