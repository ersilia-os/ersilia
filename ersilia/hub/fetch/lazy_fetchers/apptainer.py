import json
import os

from .... import ErsiliaBase, throw_ersilia_exception
from ....default import (
    API_SCHEMA_FILE,
    APIS_LIST_FILE,
    INFORMATION_FILE,
    PREDEFINED_COLUMN_FILE,
    PREDEFINED_EXAMPLE_FILES,
)
from ....setup.requirements.apptainer import ApptainerRequirement
from ....utils import apptainer
from ....utils.echo import spinner
from ....utils.exceptions_utils.cli_exceptions import ApptainerDownloadError
from .. import STATUS_FILE
from ..register.register import ModelRegisterer


class ModelApptainerFetcher(ErsiliaBase):
    """
    Fetches a model as an Apptainer image (SIF) from Ersilia's SIF bucket.

    The image is downloaded to ``SIF_DIR`` (``~/eos/sifs`` unless
    ``ERSILIA_SIF_DIR`` is set). The model's metadata is copied out of the
    image into its folder under ``~/eos/dest``, as for DockerHub models.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings.
    version : str, optional
        Image version, e.g. "v1". The newest one if not given.

    Examples
    --------
    .. code-block:: python

        fetcher = ModelApptainerFetcher(
            config_json=config
        )
        await fetcher.fetch("eos4e40")
    """

    def __init__(self, config_json=None, version=None):
        super().__init__(config_json=config_json, credentials_json=None)
        self.version = version

    def _choose_runner(self, binary, sif):
        # Plain 'apptainer exec' works on most machines. Some (e.g. Colab,
        # which runs as root) need it wrapped in 'unshare -r'.
        runner = apptainer.SimpleApptainer(binary)
        if runner.works(sif):
            return runner
        runner = apptainer.SimpleApptainer(binary, use_unshare=True)
        if runner.works(sif):
            self.logger.debug("Apptainer needs 'unshare -r' on this machine")
            return runner
        raise ApptainerDownloadError(
            os.path.basename(sif), "the image was downloaded but cannot be run"
        )

    def _copy_from_image(self, runner, sif, bundle, model_id, relative_path):
        source = runner.find_bundle_file(sif, bundle, relative_path)
        if source is None:
            self.logger.debug("{0} is not in the image".format(relative_path))
            return False
        content = runner.read_file(sif, source)
        if content is None:
            return False
        destination = os.path.join(self._model_path(model_id), relative_path)
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        with open(destination, "wb") as f:
            f.write(content)
        return True

    def _set_up(self, runner, sif, bundle, model_id):
        for relative_path in (INFORMATION_FILE, API_SCHEMA_FILE, STATUS_FILE):
            self._copy_from_image(runner, sif, bundle, model_id, relative_path)
        for relative_path in list(PREDEFINED_EXAMPLE_FILES) + [PREDEFINED_COLUMN_FILE]:
            self._copy_from_image(runner, sif, bundle, model_id, relative_path)
        self._modify_information(model_id, sif)
        # All ersilia-pack models expose one API, 'run'.
        with open(
            os.path.join(self._get_bundle_location(model_id), APIS_LIST_FILE), "w"
        ) as f:
            f.write("run" + os.linesep)

    def _modify_information(self, model_id, sif):
        information_file = os.path.join(self._model_path(model_id), INFORMATION_FILE)
        try:
            with open(information_file, "r") as f:
                data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            self.logger.error("Information file not found, not modifying anything")
            return
        data["service_class"] = "apptainer"
        data["size"] = os.path.getsize(sif) / (1024 * 1024)
        data["apptainer_version"] = self.version
        with open(information_file, "w") as f:
            json.dump(data, f, indent=4)

    @throw_ersilia_exception()
    async def fetch(self, model_id: str):
        """
        Fetch the model as an Apptainer image.

        Parameters
        ----------
        model_id : str
            ID of the model.

        Raises
        ------
        ApptainerNotLinuxError
            If this machine does not run Linux.
        ApptainerNotInstalledError
            If Apptainer is not installed.
        ApptainerImageNotFoundError
            If the model, or the requested version, has no image.
        """
        # Checked before anything is downloaded or written.
        binary = ApptainerRequirement().check()
        if self.version is None:
            self.version = apptainer.latest_version(model_id)
        verbose = getattr(self.logger, "verbosity", 0) == 1
        sif = apptainer.download(model_id, self.version, verbose=verbose)
        runner = self._choose_runner(binary, sif)
        bundle = runner.find_bundle(sif)
        if bundle is None:
            raise ApptainerDownloadError(
                model_id, "the image does not contain a model bundle"
            )
        data = {
            "version": self.version,
            "sif_path": sif,
            "bundle_path": bundle,
            "binary": binary,
            "use_unshare": runner.use_unshare,
        }
        mr = ModelRegisterer(model_id=model_id, config_json=self.config_json)
        await mr.register(apptainer=data)
        spinner(
            "Setting up the model",
            self._set_up,
            runner,
            sif,
            bundle,
            model_id,
            done="Model set up.",
        )
