import os
import platform


def is_inside_docker():
    """
    Check if the current environment is inside a Docker container.

    Returns
    -------
    bool
        True if inside a Docker container, False otherwise.
    """
    if os.path.isfile("/.dockerenv"):
        return True
    else:
        return False


class SystemChecker(object):
    """
    A class to check various system properties.

    Methods
    -------
    is_arm64()
        Check if the system architecture is ARM64.
    is_github_action()
        Check if the code is running inside a GitHub Action.
    is_inside_docker()
        Check if the code is running inside a Docker container.
    """

    def __init__(self):
        self.uname = platform.uname()

    def is_arm64(self):
        """
        Check if the system architecture is ARM64.

        Returns
        -------
        bool
            True if the system architecture is ARM64, False otherwise.
        """
        if self.uname.machine in ["arm64", "aarch64"]:
            return True
        if "arm64" in self.uname.version.lower():
            return True
        if "aarch64" in self.uname.version.lower():
            return True
        return False

    def is_github_action(self):
        """
        Check if the code is running inside a GitHub Action.

        Returns
        -------
        bool
            True if running inside a GitHub Action, False otherwise.
        """
        if os.environ.get("GITHUB_ACTIONS") == "true":
            return True
        else:
            return False

    def is_inside_docker(self):
        """
        Check if the code is running inside a Docker container.

        Returns
        -------
        bool
            True if running inside a Docker container, False otherwise.
        """
        return is_inside_docker()
