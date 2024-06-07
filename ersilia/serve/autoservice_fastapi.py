import os
import json

from .. import ErsiliaBase

from .services_fastapi import (
    SystemBundleService,
    CondaEnvironmentService,
    DockerImageService,
    PulledDockerImageService,
    HostedService,
)

from ..default import SERVICE_CLASS_FILE, IS_FETCHED_FROM_DOCKERHUB_FILE, IS_FETCHED_FROM_HOSTED_FILE, DEFAULT_BATCH_SIZE


DEFAULT_OUTPUT = None


class AutoServiceFastAPI(ErsiliaBase):
    def __init__(
        self,
        model_id,
        service_class=None,
        config_json=None,
        preferred_port=None,
        url=None,
    ):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.logger.debug("Setting FastAPI AutoService for {0}".format(model_id))
        self.config_json = config_json
        self.model_id = model_id
        self._preferred_port = preferred_port
        self._url = url
        self._service_class = service_class
        if service_class is None:
            self.logger.debug("No service class provided, deciding automatically")
            service_class_file = os.path.join(
                self._get_bundle_location(model_id), SERVICE_CLASS_FILE
            )
            s = None
            if os.path.exists(service_class_file):
                self.logger.debug(
                    "Service class file exists in folder {0}".format(service_class_file)
                )
                with open(service_class_file, "r") as f:
                    s = f.read()
                if not s:
                    s = None
            if s is not None:
                self.logger.debug("Service class: {0}".format(s))
                if s == "system":
                    self.service = SystemBundleService(
                        model_id, config_json=config_json, preferred_port=preferred_port
                    )
                elif s == "conda":
                    self.service = CondaEnvironmentService(
                        model_id, config_json=config_json, preferred_port=preferred_port
                    )
                elif s == "docker":
                    self.service = DockerImageService(
                        model_id, config_json=config_json, preferred_port=preferred_port
                    )
                elif s == "pulled_docker":
                    self.service = PulledDockerImageService(
                        model_id, config_json=config_json, preferred_port=preferred_port
                    )
                elif s == "hosted":
                    self.service = HostedService(
                        model_id, config_json=config_json, url=url
                    )
                else:
                    self.service = None
                self._service_class = s

    def _was_fetched_from_dockerhub(self):
        from_dockerhub_file = os.path.join(
            self._dest_dir, self.model_id, IS_FETCHED_FROM_DOCKERHUB_FILE
        )
        if not os.path.exists(from_dockerhub_file):
            return False
        with open(from_dockerhub_file, "r") as f:
            data = json.load(f)
            return data["docker_hub"]

    def _was_fetched_from_hosted_url(self):
        from_hosted_file = os.path.join(
            self._dest_dir, self.model_id, IS_FETCHED_FROM_HOSTED_FILE
        )
        if not os.path.exists(from_hosted_file):
            return False
        with open(from_hosted_file, "r") as f:
            data = json.load(f)
            return data["hosted_url"]

    def _set_api(self, api_name):
        def _method(input, output=DEFAULT_OUTPUT, batch_size=DEFAULT_BATCH_SIZE):
            return self.api(api_name, input, output, batch_size)

        setattr(self, api_name, _method)

    def _set_apis(self):
        if self.service is None:
            return
        self.logger.debug("Setting APIs for {0}".format(self.model_id))
        self.apis_list = self.service.get_apis()
    
    def _service_class_loader(self, service_class):
        if type(service_class) is SystemBundleService:
            self._service_class = "system"
            return service_class
        elif type(service_class) is CondaEnvironmentService:
            self._service_class = "conda"
            return service_class
        elif type(service_class) is DockerImageService:
            self._service_class = "docker"
            return service_class
        elif type(service_class) is PulledDockerImageService:
            self._service_class = "pulled_docker"
            return service_class
        elif type(service_class) is HostedService:
            self._service_class = "hosted"
            return service_class
        else:
            self._service_class = service_class
            if service_class == "system":
                return SystemBundleService
            elif service_class == "conda":
                return CondaEnvironmentService
            elif service_class == "docker":
                return DockerImageService
            elif service_class == "pulled_docker":
                return PulledDockerImageService
            elif service_class == "hosted":
                return HostedService
            raise Exception()
        
    def api(
        self, api_name, input, output=DEFAULT_OUTPUT, batch_size=DEFAULT_BATCH_SIZE
    ):
        self.logger.debug("API: {0}".format(api_name))
        self.logger.debug("MODEL ID: {0}".format(self.model_id))
        self.logger.debug("SERVICE URL: {0}".format(self.service.url))
        if batch_size is None:
            batch_size = DEFAULT_BATCH_SIZE
        else:
            batch_size = batch_size
        _api = FastAPIAppInterface(
            model_id=self.model_id,
            url=self.service.url,
            api_name=api_name,
            save_to_lake=False,
            config_json=self.config_json,
        )
        self.logger.debug("Meta: {0}".format(self._meta))
        for result in _api.post(input=input, output=output, batch_size=batch_size):
            if self._meta is None:
                do_meta = True
            else:
                if api_name not in self._meta:
                    do_meta = True
                else:
                    do_meta = False
            if do_meta:
                self.logger.debug("Metadata needs to be calculated")
                self._latest_meta = _api.meta()
                self._meta = {api_name: self._latest_meta}
            else:
                pass
            if api_name not in self._meta:
                self._meta = {api_name: _api.meta()}
            yield result
            

