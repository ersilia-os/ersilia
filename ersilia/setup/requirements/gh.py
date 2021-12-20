class GithubCliRequirement(object):
    def __ini__(self):
        self.name = "gh"
        if not self.is_installed():
            self.install()

    def is_installed(self):
        pass

    def install(self):
        pass
