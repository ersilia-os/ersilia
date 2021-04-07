from .services import SystemBundleService, CondaEnvironmentService, DockerImageService
from .. import ErsiliaBase
import os

class AutoService(ErsiliaBase):

    def __init__(self, model_id, service_class=None, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.model_id = model_id
        if service_class is None:
            # decide automatically
            service_class_file=os.path.join(self._get_bundle_location(model_id),"service_class.txt")
            if os.path.exists(service_class_file):
                with open(service_class_file, "r") as f:
                    s=f.read()
                if s == "system":
                    self.service=SystemBundleService(model_id, config_json=config_json)
                elif s == "conda":
                    self.service = CondaEnvironmentService(model_id, config_json=config_json)
                elif s == "docker":
                    self.service = DockerImageService(model_id, config_json=config_json)
                else:
                    self.service = None
            else:
                with open(service_class_file, "w") as f:
                    if SystemBundleService(model_id, config_json=config_json).is_available():
                        self.service = SystemBundleService(model_id, config_json=config_json)
                        f.write("system")
                    elif CondaEnvironmentService(model_id, config_json=config_json).is_available():
                        self.service = CondaEnvironmentService(model_id, config_json=config_json)
                        f.write("conda")
                    elif DockerImageService(model_id, config_json=config_json).is_available():
                        self.service = DockerImageService(model_id, config_json=config_json)
                        f.write("docker")
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
        apis_list= os.path.join(self._get_bundle_location(self.model_id),"apis_list.txt")
        if os.path.exists(apis_list):
            with open(apis_list, "r") as f:
                for l in f:
                    api_name = l.rstrip()
                    self._set_api(api_name)
        else:
            with open(apis_list,"w") as f:
                for api_name in self.service._get_apis_from_bento():
                    self._set_api(api_name)
                    f.write(api_name+os.linesep)
        self.apis_list=apis_list

    def get_apis(self):
        apis=[]
        with open(self.apis_list,"r") as f:
            for l in f:
                api=l.rstrip()
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

    def api(self, api_name, input):
        return self.service.api(api_name=api_name, input=input)
