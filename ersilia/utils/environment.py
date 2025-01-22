"""
Utility functions to get information about the working environment.
"""

from importlib.metadata import distributions


class Environment(object):
    """
    Class to handle environment settings.
    """

    def __init__(self):
        self.python_packages = {dist.metadata["Name"] for dist in distributions()}

    def has_module(self, module_name):
        """Check if Python module is installed."""
        if module_name in self.python_packages:
            return True
        else:
            return False
