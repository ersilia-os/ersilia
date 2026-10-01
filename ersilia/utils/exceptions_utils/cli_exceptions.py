from .exceptions import ErsiliaError

# ruff: noqa: D101, D102

# Errors for situations a user can run into with the CLI. Each one says what
# happened and, in the hint, what to do next.


class RunNotSupportedError(ErsiliaError):
    def __init__(self, model_id, reason):
        self.message = "Model {0} cannot be run with 'ersilia run': {1}.".format(
            model_id, reason
        )
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)


class ModelNotRespondingError(ErsiliaError):
    def __init__(self, model_id, url, during_run=False):
        if during_run:
            self.message = "Model {0} stopped responding during the run.".format(
                model_id
            )
            self.hints = "No output was written. Serve it again with 'ersilia serve {0}' and rerun.".format(
                model_id
            )
        else:
            self.message = "Model {0} is not responding at {1}.".format(model_id, url)
            self.hints = (
                "It may have stopped. Serve it again with 'ersilia serve {0}'.".format(
                    model_id
                )
            )
        ErsiliaError.__init__(self, self.message, self.hints)


class ResultCountMismatchError(ErsiliaError):
    def __init__(self, model_id, got, expected):
        self.message = "Model {0} returned {1} results for {2} inputs, so they cannot be matched to the inputs.".format(
            model_id, got, expected
        )
        self.hints = "No output was written. Please report this model at https://github.com/ersilia-os/ersilia/issues."
        ErsiliaError.__init__(self, self.message, self.hints)


class NoInputProcessedError(ErsiliaError):
    def __init__(self, total):
        if total == 1:
            self.message = "The input could not be processed."
        else:
            self.message = "None of the {0:,} inputs could be processed.".format(total)
        self.hints = "No output was written. Check that the inputs are valid for this model (e.g. valid SMILES)."
        ErsiliaError.__init__(self, self.message, self.hints)


class ModelStartError(ErsiliaError):
    def __init__(self, model_id, reason, hint):
        self.message = "Model {0} {1}.".format(model_id, reason)
        self.hints = hint
        ErsiliaError.__init__(self, self.message, self.hints)


class HubUnreachableError(ErsiliaError):
    def __init__(self):
        self.message = "Could not reach the Ersilia Model Hub."
        self.hints = "Check your internet connection and try again. You can also browse the models at https://catalog.ersilia.io"
        ErsiliaError.__init__(self, self.message, self.hints)


class DockerNotInstalledError(ErsiliaError):
    def __init__(self):
        self.message = "Docker is not installed."
        self.hints = "Install it from https://docs.docker.com/get-docker/ and try again, or fetch the model with --from_github."
        import platform

        if platform.system() == "Linux":
            self.hints += "\nOn Linux without Docker (e.g. HPC clusters), you can use Apptainer instead: --from_apptainer."
        ErsiliaError.__init__(self, self.message, self.hints)


class PortInUseError(ErsiliaError):
    def __init__(self, port):
        self.message = "Port {0} is already in use.".format(port)
        self.hints = "Choose another with --port, or leave it out to use a free port."
        ErsiliaError.__init__(self, self.message, self.hints)


class ImageMissingError(ErsiliaError):
    def __init__(self, model_id, image):
        self.message = "The Docker image of model {0} ({1}) is missing.".format(
            model_id, image
        )
        self.hints = "It may have been removed with 'docker rmi'. Run 'ersilia delete {0}' and then 'ersilia fetch {0}'.".format(
            model_id
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class ImagePullError(ErsiliaError):
    def __init__(self, model_id, reason, hint):
        self.message = "Could not download the Docker image of model {0}: {1}.".format(
            model_id, reason
        )
        self.hints = hint
        ErsiliaError.__init__(self, self.message, self.hints)


class ConfigFileError(ErsiliaError):
    def __init__(self, path):
        self.message = "Ersilia's settings file {0} is damaged.".format(path)
        self.hints = "Delete it and run the command again; it will be created anew."
        ErsiliaError.__init__(self, self.message, self.hints)


class ApptainerNotLinuxError(ErsiliaError):
    def __init__(self, system):
        self.message = "Apptainer images can only be run on Linux, not on {0}.".format(
            system
        )
        self.hints = "Fetch the model from DockerHub instead (the default): leave out --from_apptainer."
        ErsiliaError.__init__(self, self.message, self.hints)


APPTAINER_INSTALL_HINT = (
    "On an HPC cluster it is usually available as a module: run 'module load apptainer' "
    "(or 'module load singularity'), or ask your administrators.\n"
    "To install it yourself, see https://apptainer.org/docs/admin/main/installation.html "
    "(on Ubuntu: sudo add-apt-repository -y ppa:apptainer/ppa && sudo apt install -y apptainer)."
)


class ApptainerNotInstalledError(ErsiliaError):
    def __init__(self, model_id=None):
        if model_id:
            self.message = (
                "Model {0} runs with Apptainer, but Apptainer is not available.".format(
                    model_id
                )
            )
        else:
            self.message = "Apptainer is not installed."
        self.hints = APPTAINER_INSTALL_HINT
        if not model_id:
            self.hints += "\nOr fetch the model from DockerHub instead: leave out --from_apptainer."
        ErsiliaError.__init__(self, self.message, self.hints)


class ApptainerNotUsableError(ErsiliaError):
    def __init__(self, model_id, detail=None):
        self.message = "Apptainer is installed but cannot run model {0}'s image on this machine.".format(
            model_id
        )
        hints = []
        if detail:
            hints.append("Apptainer said: {0}".format(detail))
        hints.append(
            "It usually needs unprivileged user namespaces, or a setuid installation "
            "(apptainer-suid). On a cluster, ask your administrators; see "
            "https://apptainer.org/docs/admin/main/user_namespace.html"
        )
        hints.append(
            "The image is kept, so fetching again once Apptainer works does not download it again."
        )
        self.hints = "\n".join(hints)
        ErsiliaError.__init__(self, self.message, self.hints)


class ApptainerImageNotFoundError(ErsiliaError):
    def __init__(self, model_id, version=None):
        if version is None:
            self.message = "Model {0} has no Apptainer image.".format(model_id)
            self.hints = "Fetch it from DockerHub instead: 'ersilia fetch {0}'.".format(
                model_id
            )
        else:
            self.message = "Model {0} has no Apptainer image of version {1}.".format(
                model_id, version
            )
            self.hints = "Leave out --version to get the latest one."
        ErsiliaError.__init__(self, self.message, self.hints)


class ApptainerDownloadError(ErsiliaError):
    def __init__(self, model_id, reason):
        self.message = (
            "Could not download the Apptainer image of model {0}: {1}.".format(
                model_id, reason
            )
        )
        self.hints = "Check your internet connection and disk space, and try again."
        ErsiliaError.__init__(self, self.message, self.hints)
