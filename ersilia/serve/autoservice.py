from .services import SystemBundleService, CondaEnvironmentService, DockerImageService


class AutoService(object):

    def __init__(self, model_id, service_class=None, config_json=None):
        self.model_id = model_id
        if service_class is None:
            # decide automatically
            if SystemBundleService(model_id, config_json=config_json).is_available():
                self.service = SystemBundleService(model_id, config_json=config_json)
            elif CondaEnvironmentService(model_id, config_json=config_json).is_available():
                self.service = CondaEnvironmentService(model_id, config_json=config_json)
            elif DockerImageService(model_id, config_json=config_json).is_available():
                self.service = DockerImageService(model_id, config_json=config_json)
            else:
                self.service = None
        else:
            # predefined service class
            if service_class(model_id, config_json).is_available():
                self.service = service_class(model_id, config_json=config_json)
            else:
                self.service = None
        self._set_apis()

    def _set_api(self, api_name):
        def _method(x):
            return self.service.api(api_name, x)
        setattr(self, api_name, _method)

    def _set_apis(self):
        if self.service is None:
            return
        for api_name in self.service._get_apis_from_bento():
            self._set_api(api_name)

    def is_available(self):
        if self.service is None:
            return False
        else:
            return True

    def serve(self):
        self.service.serve()

    def close(self):
        self.service.close()

    def api(self, api_name, input):
        return self.service.api(api_name=api_name, input=input)
