from .exceptions import ErsiliaError


class ErsiliaHomePathError(ErsiliaError):
    def __init__(self):
        pass


class GitLfsSetupError(ErsiliaError):
    def __init__(
        self,
    ):
        self.message = self._get_message()
        self.hints = self._get_hints()
        super().__init__(self.message, self.hints)

    def _get_message(self):
        text = "Git LFS is not installed! Git LFS is needed to download large files from our GitHub Model Hub repositories."
        return text

    def _get_hints(self):
        text = "An easy way to install Git LFS is the following Conda command:\n"
        text += "$ conda install -c conda-forge git-lfs"
        return text
