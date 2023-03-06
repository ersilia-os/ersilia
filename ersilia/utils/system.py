import platform


class SystemChecker(object):
    def __init__(self):
        self.uname = platform.uname()

    def is_arm64(self):
        if self.uname.machine == "arm64":
            return True
        if "arm64" in self.uname.version.lower():
            return True
        return False
