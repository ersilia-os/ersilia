import asyncio
import json
import os

from .... import EOS, ErsiliaBase, throw_ersilia_exception
from ....default import (
    API_SCHEMA_FILE,
    DOCKERHUB_LATEST_TAG,
    DOCKERHUB_ORG,
    INFORMATION_FILE,
    PACK_METHOD_BENTOML,
    PACK_METHOD_FASTAPI,
    PREDEFINED_EXAMPLE_FILES,
)
from ....serve.services import PulledDockerImageService
from ....setup.requirements.docker import DockerRequirement
from ....utils.docker import SimpleDocker
from ....utils.resolvers import PackMethodResolver
from ...pull.pull import ModelPuller
from .. import STATUS_FILE
from ..actions.template_resolver import TemplateResolver
from ..register.register import ModelRegisterer


class ModelDockerHubFetcher(ErsiliaBase):
    """
    A class used to fetch models from DockerHub.

    Attributes
    ----------
    overwrite : bool
        Flag to indicate whether to overwrite existing files.
    config_json : dict
        Configuration settings in JSON format.
    simple_docker : SimpleDocker
        Instance of SimpleDocker for Docker operations.

    Methods
    -------
    is_docker_installed()
        Check if Docker is installed.
    is_docker_active()
        Check if Docker is active.
    is_available(model_id)
        Check if the model is available locally or on DockerHub.
    write_apis(model_id)
        Write APIs for the model.
    copy_information(model_id)
        Copy the information file from the model container.
    copy_metadata(model_id)
        Copy the metadata file from the model container.
    copy_status(model_id)
        Copy the status file from the model container.
    copy_example_if_available(model_id)
        Copy example files from the model container if available.
    modify_information(model_id)
        Modify the information file copied from the model container.
    fetch(model_id)
        Fetch the model from DockerHub.
    """

    def __init__(
        self,
        overwrite=None,
        config_json=None,
        img_tag=None,
        force_with_bentoml=False,
        force_with_fastapi=False,
    ):
        super().__init__(config_json=config_json, credentials_json=None)
        self.simple_docker = SimpleDocker()
        self.overwrite = overwrite
        self.img_tag = img_tag or DOCKERHUB_LATEST_TAG
        self.pack_method = None
        self.force_with_bentoml = force_with_bentoml
        self.force_with_fastapi = force_with_fastapi

    async def _resolve_pack_method(self, model_id):
        if self.force_with_fastapi:
            pack_method = PACK_METHOD_FASTAPI
        elif self.force_with_bentoml:
            pack_method = PACK_METHOD_BENTOML
        else:
            self.logger.debug("Resolving pack method from GitHub repository...")
            pmr = PackMethodResolver(model_id=model_id)
            pack_method = pmr.resolve_pack_method_from_github_metadata()
            if not pack_method:
                self.logger.debug("Resolving pack method with TemplateResolver...")
                tr = TemplateResolver(model_id=model_id)
                if tr.is_fastapi():
                    pack_method = PACK_METHOD_FASTAPI
                elif tr.is_bentoml():
                    pack_method = PACK_METHOD_BENTOML
                else:
                    raise Exception("Pack method could not be resolved")
        self.logger.debug(f"Pack method resolved: {pack_method}")
        self.pack_method = pack_method

    def is_docker_installed(self) -> bool:
        """
        Check if Docker is installed.

        Returns
        -------
        bool
            True if Docker is installed, False otherwise.
        """
        return DockerRequirement().is_installed()

    def is_docker_active(self) -> bool:
        """
        Check if Docker is active.

        Returns
        -------
        bool
            True if Docker is active, False otherwise.
        """
        return DockerRequirement().is_active()

    def is_available(self, model_id: str) -> bool:
        """
        Check if the model is available locally or on DockerHub.

        Parameters
        ----------
        model_id : str
            ID of the model to check.

        Returns
        -------
        bool
            True if the model is available, False otherwise.
        """
        mp = ModelPuller(
            model_id=model_id,
            overwrite=self.overwrite,
            config_json=self.config_json,
            docker_tag=self.img_tag,
        )
        if mp.is_available_locally():
            return True
        if mp.is_available_in_dockerhub():
            return True
        return False

    async def write_apis(self, model_id: str):
        """
        Write APIs for the model.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        self.logger.debug("Writing APIs")
        di = PulledDockerImageService(
            model_id=model_id, config_json=self.config_json, preferred_port=None
        )
        di.serve()
        di.close()

    async def _copy_from_bentoml_image(self, model_id: str, file: str):
        """
        Copy a file from a BentoML image.

        Parameters
        ----------
        model_id : str
            ID of the model.
        file : str
            Name of the file to copy.
        """
        fr_file = f"/root/eos/dest/{model_id}/{file}"
        to_file = f"{EOS}/dest/{model_id}/{file}"
        try:
            await self.simple_docker.cp_from_image(
                img_path=fr_file,
                local_path=to_file,
                org=DOCKERHUB_ORG,
                img=model_id,
                tag=self.img_tag,
            )
        except Exception as e:
            self.logger.error(f"Exception when copying: {e}")

    async def _copy_from_ersiliapack_image(self, model_id: str, file: str):
        """
        Copy a file from an ErsiliaPack image.

        Parameters
        ----------
        model_id : str
            ID of the model.
        file : str
            Name of the file to copy.
        """
        fr_file = f"/root/{file}"
        to_file = f"{EOS}/dest/{model_id}/{file}"
        await self.simple_docker.cp_from_image(
            img_path=fr_file,
            local_path=to_file,
            org=DOCKERHUB_ORG,
            img=model_id,
            tag=self.img_tag,
        )

    async def _copy_from_image_to_local(self, model_id: str, file: str):
        """
        Copy a file from the Docker image to the local filesystem.

        Parameters
        ----------
        model_id : str
            ID of the model.
        file : str
            Name of the file to copy.
        """
        if not self.pack_method:
            raise Exception("Pack method not resolved")
        if self.pack_method == PACK_METHOD_BENTOML:
            await self._copy_from_bentoml_image(model_id, file)
        elif self.pack_method == PACK_METHOD_FASTAPI:
            await self._copy_from_ersiliapack_image(model_id, file)
        else:
            raise Exception(
                f"Unresolved method, should be either {PACK_METHOD_FASTAPI} or {PACK_METHOD_BENTOML}"
            )

    async def copy_information(self, model_id: str):
        """
        Copy the information file from the model container.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        self.logger.debug("Copying information file from model container")
        await self._copy_from_image_to_local(model_id, INFORMATION_FILE)

    async def copy_metadata(self, model_id: str):
        """
        Copy the metadata file from the model container.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        self.logger.debug("Copying api_schema_file file from model container")
        await self._copy_from_image_to_local(model_id, API_SCHEMA_FILE)

    async def copy_status(self, model_id: str):
        """
        Copy the status file from the model container.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        self.logger.debug("Copying status file from model container")
        await self._copy_from_image_to_local(model_id, STATUS_FILE)

    async def copy_example_if_available(self, model_id: str):
        """
        Copy example files from the model container if available.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        for pf in PREDEFINED_EXAMPLE_FILES:
            await self._copy_from_image_to_local(model_id, pf)

    async def modify_information(self, model_id: str):
        """
        Modify the information file copied from the model directory.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        information_file = os.path.join(self._model_path(model_id), INFORMATION_FILE)
        mp = ModelPuller(
            model_id=model_id, config_json=self.config_json, docker_tag=self.img_tag
        )
        try:
            with open(information_file, "r") as infile:
                data = json.load(infile)
        except FileNotFoundError:
            self.logger.error("Information file not found, not modifying anything")
            return None

        data["service_class"] = "pulled_docker"
        data["size"] = mp._get_size_of_local_docker_image_in_mb()
        with open(information_file, "w") as outfile:
            json.dump(data, outfile, indent=4)

    @throw_ersilia_exception()
    async def fetch(self, model_id: str):
        """
        Fetch the model from DockerHub.

        Parameters
        ----------
        model_id : str
            ID of the model.
        """
        await self._resolve_pack_method(model_id)
        mp = ModelPuller(
            model_id=model_id, config_json=self.config_json, docker_tag=self.img_tag
        )
        self.logger.debug("Pulling model image from DockerHub")
        await mp.async_pull()
        mr = ModelRegisterer(model_id=model_id, config_json=self.config_json)
        self.logger.debug("Asynchronous and concurrent execution started!")
        await asyncio.gather(
            mr.register(is_from_dockerhub=True, img_tag=self.img_tag),
            self.write_apis(model_id),
            self.copy_information(model_id),
            self.modify_information(model_id),
            self.copy_metadata(model_id),
            self.copy_status(model_id),
            self.copy_example_if_available(model_id),
        )
