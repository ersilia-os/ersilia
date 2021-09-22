import os

from .services import (
    SystemBundleService,
    VenvEnvironmentService,
    CondaEnvironmentService,
    DockerImageService,
)
from .api import Api
from ..default import DEFAULT_BATCH_SIZE
from .. import ErsiliaBase

DEFAULT_OUTPUT = None
DEFAULT_BATCH_SIZE = None


class AutoService(ErsiliaBase):
    def __init__(self, model_id, service_class=None, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.logger.debug("Setting AutoService for {0}".format(model_id))
        self.config_json = config_json
        self.model_id = model_id
        if service_class is None:
            self.logger.debug("No service class provided, deciding automatically")
            # decide automatically
            service_class_file = os.path.join(
                self._get_bundle_location(model_id), "service_class.txt"
            )
            if os.path.exists(service_class_file):
                self.logger.debug(
                    "Service class file exists in folder {0}".format(service_class_file)
                )
                with open(service_class_file, "r") as f:
                    s = f.read()
                self.logger.debug("Service class: {0}".format(s))
                if s == "system":
                    self.service = SystemBundleService(
                        model_id, config_json=config_json
                    )
                elif s == "venv":
                    self.service = VenvEnvironmentService(
                        model_id, config_json=config_json
                    )
                elif s == "conda":
                    self.service = CondaEnvironmentService(
                        model_id, config_json=config_json
                    )
                elif s == "docker":
                    self.service = DockerImageService(model_id, config_json=config_json)
                else:
                    self.service = None
            else:
                self.logger.debug(
                    "No service class file exists in {0}".format(service_class_file)
                )
                with open(service_class_file, "w") as f:
                    if SystemBundleService(
                        model_id, config_json=config_json
                    ).is_available():
                        self.service = SystemBundleService(
                            model_id, config_json=config_json
                        )
                        f.write("system")
                    elif VenvEnvironmentService(
                        model_id, config_json=config_json
                    ).is_available():
                        self.service = VenvEnvironmentService(
                            model_id, config_json=config_json
                        )
                        f.write("venv")
                    elif CondaEnvironmentService(
                        model_id, config_json=config_json
                    ).is_available():
                        self.service = CondaEnvironmentService(
                            model_id, config_json=config_json
                        )
                        f.write("conda")
                    elif DockerImageService(
                        model_id, config_json=config_json
                    ).is_available():
                        self.service = DockerImageService(
                            model_id, config_json=config_json
                        )
                        f.write("docker")
                    else:
                        self.service = None
        else:
            self.logger.info("Service class provided")
            # predefined service class
            if service_class(model_id, config_json).is_available():
                self.service = service_class(model_id, config_json=config_json)
            else:
                self.service = None
        self._set_apis()

    def _set_api(self, api_name):
        def _method(input, output=DEFAULT_OUTPUT, batch_size=DEFAULT_BATCH_SIZE):
            return self.api(api_name, input, output, batch_size)

        setattr(self, api_name, _method)

    def _set_apis(self):
        if self.service is None:
            return
        apis_list = os.path.join(
            self._get_bundle_location(self.model_id), "apis_list.txt"
        )
        if os.path.exists(apis_list):
            with open(apis_list, "r") as f:
                for l in f:
                    api_name = l.rstrip()
                    self._set_api(api_name)
        else:
            with open(apis_list, "w") as f:
                for api_name in self.service._get_apis_from_bento():
                    self._set_api(api_name)
                    f.write(api_name + os.linesep)
        self.apis_list = apis_list

    def get_apis(self):
        apis = []
        with open(self.apis_list, "r") as f:
            for l in f:
                api = l.rstrip()
                apis.append(api)
        return sorted(apis)

    def is_available(self):
        if self.service is None:
            return False
        else:
            return True

    def serve(self):
        self.service.serve()

    def close(self):
        self.service.close()

    def api(self, api_name, input, output=DEFAULT_OUTPUT, batch_size=DEFAULT_BATCH_SIZE):
        self.logger.debug("API: {0}".format(api_name))
        self.logger.debug("MODEL ID: {0}".format(self.model_id))
        self.logger.debug("SERVICE URL: {0}".format(self.service.url))
        if batch_size is None:
            batch_size = DEFAULT_BATCH_SIZE
        else:
            batch_size = batch_size
        _api = Api(
            model_id=self.model_id,
            url=self.service.url,
            api_name=api_name,
            config_json=self.config_json,
        )
        for result in _api.post(input=input, output=output, batch_size=batch_size):
            yield result
