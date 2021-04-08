"""
Utility functions to get information about the working environment.
"""
import pkg_resources


class Environment(object):
    def __init__(self):
        self.python_packages = {pkg.key for pkg in pkg_resources.working_set}

    def has_module(self, module_name):
        """Check if Python module is installed."""
        if module_name in self.python_packages:
            return True
        else:
            return False
