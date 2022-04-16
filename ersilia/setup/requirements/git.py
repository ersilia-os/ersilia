class GithubCliRequirement(object):
    def __ini__(self):
        self.name = "gh"
        if not self.is_installed():
            self.install()

    def is_installed(self):
        pass

    def install(self):
        pass


class GitLfsRequirement(object):
    def __init__(self):
        pass

    def is_installed(self):
        pass