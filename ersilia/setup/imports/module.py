"""
from pkgutil import iter_modules


class Module(object):

    def __init__(self, BaseModule):
        self.module = BaseModule()

    def exists(self):
        return self.module.name in (name for loader, name, ispkg in iter_modules())

    def install(self):
        self.module.install()

    def uninstall(self):
        self.module.uninstall()

    def import(self):
        if not self.exists(self.module.name):
            pass
        else:
            pass
"""
