import sys
import subprocess as s


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
        try:
            check = s.run(["git-lfs"], capture_output=True)

        except ModuleNotFoundError:
            sys.exit(
                "\nGit LFS is not installed. We recommend installing Git LFS to"
                " use large size models.\n\nCheck out https://git-lfs.github.com/"
                " on how to install. After installation, simply use the command"
                " `git lfs install` in your repository.\n"
            )
