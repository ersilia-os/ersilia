from .exceptions import ErsiliaError

# ruff: noqa: D101, D102


class FetchErsiliaError(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = "Error occured while fetching model: {0}".format(self.model_id)
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)


class InvalidUrlError(ErsiliaError):
    def __init__(self, url):
        self.message = "Provided URL is invalid: {0}".format(url)
        self.hints = "Open a browser and check that the URL is valid. You should see an API interface."
        ErsiliaError.__init__(self, self.message, self.hints)


class S3DownloaderError(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = "Error occured while fetching model from S3: {0}".format(
            self.model_id
        )
        self.hints = (
            "Make sure that this model is actually in Ersilia Model Hub's S3 bucket"
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class GetFetchErsiliaError(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = "Error occured while fetching model: {0}".format(self.model_id)
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)


class FolderNotFoundError(ErsiliaError):
    def __init__(self, folder_name):
        self.folder_name = folder_name
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Folder " + self.folder_name + " not found!"
        return text

    def _get_hints(self):
        text = "Check that you can clone repositories from github\n"
        return text


class CondaEnvironmentExistsError(ErsiliaError):
    def __init__(self, env_name):
        self.environment_name = env_name
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Environment " + self.environment_name + " does not exist!"
        return text

    def _get_hints(self):
        text = "Try to look at the log files and see if there is an error in the conda installation\n"
        text += (
            "Please report the error at:\n - https://github.com/ersilia-os/ersilia\n"
        )
        return text


class ModelPackageInstallError(ErsiliaError):
    def __init__(self, package_name):
        self.package_name = package_name
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = (
            'Error occured while installing package by running "'
            + self.package_name
            + '" command \n'
        )
        return text

    def _get_hints(self):
        text = "Try to manually activate the model environment and install package manually.\n"
        text += "If this does not work, please report the error at:\n - https://github.com/ersilia-os/ersilia\n\n"
        return text


class VirtualEnvironmentSetupError(ErsiliaError):
    def __init__(self, venv_name):
        self.virtual_env_name = venv_name
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Virtual Environment " + self.virtual_env_name + " does not exist!"
        return text

    def _get_hints(self):
        text = "Try to look at the log files and see if there is an error in the venv installation\n"
        text += (
            "Please report the error at:\n - https://github.com/ersilia-os/ersilia\n"
        )
        return text


class OutputDataTypesNotConsistentError(ErsiliaError):
    def __init__(self):
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Output data types are not consistent"
        return text

    def _get_hints(self):
        text = "This message is related to a bad development of the model. As an end user, there is not much you can do about it. Please reach out to Ersilia directly to report this error."
        return text


class StandardModelExampleError(ErsiliaError):
    def __init__(self, model_id, file_name):
        self.model_id = model_id
        self.file_name = file_name
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Standard model run from CSV was not possible for model {0}".format(
            self.model_id
        )
        text += "\n"
        text += "Output file {0} was not created successfully".format(self.file_name)
        return text

    def _get_hints(self):
        text = "If you fetch this model from Docker Hub, or you are running it through URL, this is the first time run is executed in your local computer. Reach out to Ersilia to get specific help."
        return text


class DockerNotActiveError(ErsiliaError):
    def __init__(self):
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Cannot fetch model from Docker Hub since Docker is not active."
        return text

    def _get_hints(self):
        text = "Make sure that Docker is running on your computer. We recommend to use Docker Desktop."
        return text


class NotInstallableError(ErsiliaError):
    def __init__(self, model_id, packing_strategy):
        self.packing_strategy = packing_strategy
        self.model_id = model_id
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = f"Model {self.model_id} is not installable with {self.packing_strategy}"
        return text

    def _get_hints(self):
        text = f"This model is not compatible with {self.packing_strategy}. Please check the model structure or reach out to Ersilia directly to report this error."
        return text


class NotInstallableWithFastAPI(NotInstallableError):
    def __init__(self, model_id):
        super.__init__(model_id, "FastAPI")


class NotInstallableWithBentoML(NotInstallableError):
    def __init__(self, model_id):
        super.__init__(model_id, "BentoML")


class SniffFastApiColumnsDontMatch(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Column names in columns file and example output files are not the same! Please revise the model."
        return text

    def _get_hints(self):
        text = "Check the model repository and inspect the files to make sure they are properly constructed."
        return text


class SniffFastApiColumnTypesIncompatibility(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Column types are not purely string or purely numeric. In Ersilia Pack (FastAPI) deployments, they need to be pure."
        return text

    def _get_hints(self):
        text = "Check the model repository and make sure that the type column is purely numeric or string."
        return text


class WithToolFetchingNotWorking(ErsiliaError):
    def __init__(self, tool):
        assert tool in ["bentoml", "fastapi"]
        self.tool = tool
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Fetching with {0} did not work".format(self.tool)
        return text

    def _get_hints(self):
        text = "Check the model repository structure and make sure that all files and Python versions are correct to fetch with {0}".format(
            self.tool
        )
        return text
