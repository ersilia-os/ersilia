import os
import json
import asyncio
from ..register.register import ModelRegisterer

from .... import ErsiliaBase, throw_ersilia_exception
from .... import EOS
from ....default import (
    DOCKERHUB_ORG,
    DOCKERHUB_LATEST_TAG,
    PREDEFINED_EXAMPLE_FILES,
    INFORMATION_FILE,
    API_SCHEMA_FILE,
)

from ...pull.pull import ModelPuller
from ....serve.services import PulledDockerImageService
from ....setup.requirements.docker import DockerRequirement
from ....utils.docker import SimpleDocker, resolve_pack_method_docker, PACK_METHOD_BENTOML
from .. import STATUS_FILE

class ModelDockerHubFetcher(ErsiliaBase):
    def __init__(self, overwrite=None, config_json=None):
        super().__init__(config_json=config_json, credentials_json=None)
        self.simple_docker = SimpleDocker()
        self.overwrite = overwrite

    def is_docker_installed(self):
        return DockerRequirement().is_installed()

    def is_docker_active(self):
        return DockerRequirement().is_active()

    def is_available(self, model_id):
        mp = ModelPuller(
            model_id=model_id, overwrite=self.overwrite, config_json=self.config_json
        )
        if mp.is_available_locally():
            return True
        if mp.is_available_in_dockerhub():
            return True
        return False

    async def write_apis(self, model_id):
        self.logger.debug("Writing APIs")
        di = PulledDockerImageService(
            model_id=model_id, config_json=self.config_json, preferred_port=None
        )
        di.serve()
        di.close()

    async def _copy_from_bentoml_image(self, model_id, file):
        fr_file = f"/root/eos/dest/{model_id}/{file}"
        to_file = f"{EOS}/dest/{model_id}/{file}"
        try:
            await self.simple_docker.cp_from_image(
                img_path=fr_file,
                local_path=to_file,
                org=DOCKERHUB_ORG,
                img=model_id,
                tag=DOCKERHUB_LATEST_TAG,
            )
        except Exception as e:
            self.logger.error(f"Exception when copying: {e}")

    async def _copy_from_ersiliapack_image(self, model_id, file):
        fr_file = f"/root/{file}"
        to_file = f"{EOS}/dest/{model_id}/{file}"
        await self.simple_docker.cp_from_image(
            img_path=fr_file,
            local_path=to_file,
            org=DOCKERHUB_ORG,
            img=model_id,
            tag=DOCKERHUB_LATEST_TAG,
        )

    async def _copy_from_image_to_local(self, model_id, file):
        pack_method = resolve_pack_method_docker(model_id)
        if pack_method == PACK_METHOD_BENTOML:
            await self._copy_from_bentoml_image(model_id, file)
        else:
            await self._copy_from_ersiliapack_image(model_id, file)

    async def copy_information(self, model_id):
        self.logger.debug("Copying information file from model container")
        await self._copy_from_image_to_local(model_id, INFORMATION_FILE)

    async def copy_metadata(self, model_id):
        self.logger.debug("Copying api_schema_file file from model container")
        await self._copy_from_image_to_local(model_id, API_SCHEMA_FILE)

    async def copy_status(self, model_id):
        self.logger.debug("Copying status file from model container")
        await self._copy_from_image_to_local(model_id, STATUS_FILE)

    async def copy_example_if_available(self, model_id):
        # This needs to accommodate ersilia pack
        for pf in PREDEFINED_EXAMPLE_FILES:
            fr_file = f"/root/eos/dest/{model_id}/{pf}"
            to_file = f"{EOS}/dest/{model_id}/input.csv"
            try:
                self.simple_docker.cp_from_image(
                    img_path=fr_file,
                    local_path=to_file,
                    org=DOCKERHUB_ORG,
                    img=model_id,
                    tag=DOCKERHUB_LATEST_TAG,
                )
                return
            except:
                self.logger.debug("Could not find example file in docker image")

    async def modify_information(self, model_id):
        """
        Modify the information file being copied from docker container to the host machine.
        :param file: The model information file being copied.
        :param service_class_file: File containing the model service class.
        :size_file: File containing the size of the pulled docker image.
        """
        information_file = os.path.join(self._model_path(model_id), INFORMATION_FILE)
        mp = ModelPuller(model_id=model_id, config_json=self.config_json)
        try:
            with open(information_file, "r") as infile:
                data = json.load(infile)
        except FileNotFoundError:
            self.logger.error("Information file not found, not modifying anything")
            return None

        # Using this literal here to prevent a file read         
        # from service class file for a model fetched through DockerHub
        # since we already know the service class.
        data["service_class"] = "pulled_docker"
        data["size"] = mp._get_size_of_local_docker_image_in_mb()  # TODO this should probably be a util function 
        with open(information_file, "w") as outfile:
            json.dump(data, outfile, indent=4)

    @throw_ersilia_exception()
    async def fetch(self, model_id):
        mp = ModelPuller(model_id=model_id, config_json=self.config_json)
        self.logger.debug("Pulling model image from DockerHub")
        # Asynchronous pulling
        await mp.async_pull()
        mr = ModelRegisterer(model_id=model_id, config_json=self.config_json)
        # Asynchronous and concurent execution
        self.logger.debug("Asynchronous and concurrent execution started!")
        await asyncio.gather(
            mr.register(is_from_dockerhub=True),
            self.write_apis(model_id),
            self.copy_information(model_id),
            self.modify_information(model_id),
            self.copy_metadata(model_id),
            self.copy_status(model_id),
            self.copy_example_if_available(model_id)
        )