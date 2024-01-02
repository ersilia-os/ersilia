import platform
import os


class SystemChecker(object):
    def __init__(self):
        self.uname = platform.uname()

    def is_arm64(self):
        if self.uname.machine == "arm64":
            return True
        if "arm64" in self.uname.version.lower():
            return True
        return False

    def is_github_action(self):
        if os.environ.get("GITHUB_ACTIONS") == "true":
            return True
        else:
            return False
