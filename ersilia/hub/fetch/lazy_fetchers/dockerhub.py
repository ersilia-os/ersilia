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
from ....utils.docker import SimpleDocker
from ....utils.exceptions_utils.fetch_exceptions import DockerNotActiveError
from .. import STATUS_FILE


class ModelDockerHubFetcher(ErsiliaBase):
    def __init__(self, overwrite=None, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
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

    def write_apis(self, model_id):
        self.logger.debug("Writing APIs")
        di = PulledDockerImageService(
            model_id=model_id, config_json=self.config_json, preferred_port=None
        )
        di.serve()
        di.close()

    def copy_information(self, model_id):
        fr_file = "/root/eos/dest/{0}/{1}".format(model_id, INFORMATION_FILE)
        to_file = "{0}/dest/{1}/{2}".format(EOS, model_id, INFORMATION_FILE)
        self.simple_docker.cp_from_image(
            img_path=fr_file,
            local_path=to_file,
            org=DOCKERHUB_ORG,
            img=model_id,
            tag=DOCKERHUB_LATEST_TAG,
        )

    def copy_metadata(self, model_id):
        fr_file = "/root/eos/dest/{0}/{1}".format(model_id, API_SCHEMA_FILE)
        to_file = "{0}/dest/{1}/{2}".format(EOS, model_id, API_SCHEMA_FILE)
        self.simple_docker.cp_from_image(
            img_path=fr_file,
            local_path=to_file,
            org=DOCKERHUB_ORG,
            img=model_id,
            tag=DOCKERHUB_LATEST_TAG,
        )

    def copy_status(self, model_id):
        fr_file = "/root/eos/dest/{0}/{1}".format(model_id, STATUS_FILE)
        to_file = "{0}/dest/{1}/{2}".format(EOS, model_id, STATUS_FILE)
        self.simple_docker.cp_from_image(
            img_path=fr_file,
            local_path=to_file,
            org=DOCKERHUB_ORG,
            img=model_id,
            tag=DOCKERHUB_LATEST_TAG,
        )

    def copy_example_if_available(self, model_id):
        for pf in PREDEFINED_EXAMPLE_FILES:
            fr_file = "/root/eos/dest/{0}/{1}".format(model_id, pf)
            to_file = "{0}/dest/{1}/{2}".format(EOS, model_id, "input.csv")
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

    @throw_ersilia_exception
    def fetch(self, model_id):
        if not DockerRequirement().is_active():
            raise DockerNotActiveError()
        mp = ModelPuller(model_id=model_id, config_json=self.config_json)
        mp.pull()
        mr = ModelRegisterer(model_id=model_id, config_json=self.config_json)
        mr.register(is_from_dockerhub=True)
        self.write_apis(model_id)
        self.copy_information(model_id)
        self.copy_metadata(model_id)
        self.copy_status(model_id)
        self.copy_example_if_available(model_id)
