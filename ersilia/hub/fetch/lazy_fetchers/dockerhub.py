from ..register.register import ModelRegisterer

from .... import ErsiliaBase
from .... import EOS
from ....default import DOCKERHUB_ORG, DOCKERHUB_LATEST_TAG

from ...pull.pull import ModelPuller
from ....serve.services import PulledDockerImageService
from ....setup.requirements.docker import DockerRequirement
from ....utils.docker import SimpleDocker
from .. import STATUS_FILE


class ModelDockerHubFetcher(ErsiliaBase):
    def __init__(self, overwrite=None, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.simple_docker = SimpleDocker()
        self.overwrite = overwrite

    def is_docker_installed(self):
        return DockerRequirement().is_installed()

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
        fr_file = "/root/eos/dest/{0}/information.json".format(model_id)
        to_file = "{0}/dest/{1}/information.json".format(EOS, model_id)
        self.simple_docker.cp_from_image(
            img_path=fr_file,
            local_path=to_file,
            org=DOCKERHUB_ORG,
            img=model_id,
            tag=DOCKERHUB_LATEST_TAG,
        )

    def copy_metadata(self, model_id):
        fr_file = "/root/eos/dest/{0}/api_schema.json".format(model_id)
        to_file = "{0}/dest/{1}/api_schema.json".format(EOS, model_id)
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

    def fetch(self, model_id):
        mp = ModelPuller(model_id=model_id, config_json=self.config_json)
        mp.pull()
        mr = ModelRegisterer(model_id=model_id, config_json=self.config_json)
        mr.register(is_from_dockerhub=True)
        self.write_apis(model_id)
        self.copy_information(model_id)
        self.copy_metadata(model_id)
        self.copy_status(model_id)
