import os
import tempfile
import shutil
import json

from .services import (
    SystemBundleService,
    VenvEnvironmentService,
    CondaEnvironmentService,
    DockerImageService,
    DummyService,
    PulledDockerImageService,
    HostedService,
)
from .api import Api
from ..db.environments.managers import DockerManager
from .. import ErsiliaBase
from ..utils import tmp_pid_file

from ..default import (
    DEFAULT_BATCH_SIZE,
    SERVICE_CLASS_FILE,
    APIS_LIST_FILE,
    IS_FETCHED_FROM_DOCKERHUB_FILE,
    IS_FETCHED_FROM_HOSTED_FILE,
)

DEFAULT_OUTPUT = None


class ServiceFactory:
    SERVICE_CLASSES = {
        "system": SystemBundleService,
        "venv": VenvEnvironmentService,
        "conda": CondaEnvironmentService,
        "docker": DockerImageService,
        "pulled_docker": PulledDockerImageService,
        "hosted": HostedService,
        "dummy": DummyService,
    }

    @staticmethod
    def create_service(
        service_type, model_id, config_json, preferred_port=None, url=None
    ):
        service_class = ServiceFactory.SERVICE_CLASSES.get(service_type)
        if (
            service_class
            and service_class(model_id, config_json, preferred_port, url).is_available()
        ):
            return service_class(model_id, config_json, preferred_port, url)
        return None

    @staticmethod
    def auto_detect_service(model_id, config_json, preferred_port=None, url=None):
        for service_type, service_class in ServiceFactory.SERVICE_CLASSES.items():
            service_instance = service_class(model_id, config_json, preferred_port, url)
            if service_instance.is_available():
                return service_instance, service_type
        return DummyService(model_id, config_json, preferred_port), "dummy"


class AutoService(ErsiliaBase):
    def __init__(
        self,
        model_id,
        service_class=None,
        config_json=None,
        preferred_port=None,
        url=None,
    ):
        super().__init__(config_json=config_json)
        self.logger.debug(f"Setting BentoML AutoService for {model_id}")

        self.config_json = config_json
        self.model_id = model_id
        self._meta = None
        self._preferred_port = preferred_port
        self._url = url
        self._service_class = None

        if service_class:
            self.logger.info("Service class provided")
            loaded_service_class = self._service_class_loader(service_class)
            self.service = self._instantiate_service(loaded_service_class)
        else:
            self.service, self._service_class = self._load_or_auto_detect_service(
                model_id
            )

        self._set_apis()

    def _load_or_auto_detect_service(self, model_id):
        service_class_file = os.path.join(
            self._get_bundle_location(model_id), SERVICE_CLASS_FILE
        )

        if os.path.exists(service_class_file):
            with open(service_class_file, "r") as f:
                service_type = f.read().strip()
                self.logger.debug(f"Service class found in file: {service_type}")
                return (
                    ServiceFactory.create_service(
                        service_type,
                        self.model_id,
                        self.config_json,
                        self._preferred_port,
                        self._url,
                    ),
                    service_type,
                )

        service_instance, service_type = ServiceFactory.auto_detect_service(
            self.model_id, self.config_json, self._preferred_port, self._url
        )
        self._write_service_class_file(service_type, service_class_file)
        return service_instance, service_type

    def _instantiate_service(self, service_class):
        service_instance = service_class(
            self.model_id, self.config_json, self._preferred_port, self._url
        )
        if service_instance.is_available():
            return service_instance
        return None

    def _write_service_class_file(self, service_type, file_path):
        with open(file_path, "w") as f:
            f.write(service_type)
            self.logger.debug(f"Service class {service_type} written to file")


    def _service_class_loader(self, service_class):
        service_mapping = {
            SystemBundleService: ("system", SystemBundleService),
            VenvEnvironmentService: ("venv", VenvEnvironmentService),
            CondaEnvironmentService: ("conda", CondaEnvironmentService),
            DockerImageService: ("docker", DockerImageService),
            PulledDockerImageService: ("pulled_docker", PulledDockerImageService),
            HostedService: ("hosted", HostedService),
            "system": SystemBundleService,
            "venv": VenvEnvironmentService,
            "conda": CondaEnvironmentService,
            "docker": DockerImageService,
            "pulled_docker": PulledDockerImageService,
            "hosted": HostedService,
        }

        if isinstance(service_class, type):
            if service_class in service_mapping:
                self._service_class = service_mapping[service_class][0]
                return service_mapping[service_class][1]
        elif isinstance(service_class, str) and service_class in service_mapping:
            self._service_class = service_class
            return service_mapping[service_class]

        raise ValueError(f"Unknown service class: {service_class}")

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
        apis_list = os.path.join(
            self._get_bundle_location(self.model_id), APIS_LIST_FILE
        )
        if os.path.exists(apis_list):
            with open(apis_list, "r") as f:
                for l in f:
                    api_name = l.rstrip()
                    self._set_api(api_name)
        else:
            with open(apis_list, "w") as f:
                for api_name in self.service._get_apis_from_where_available():
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

    def is_served(self):
        tmp_file = tmp_pid_file(self.model_id)
        if os.path.exists(tmp_file):
            return True
        else:
            return False

    def _pids_from_file(self, fn):
        pids = []
        with open(fn, "r") as f:
            for l in f:
                pids += [int(l.split(" ")[0])]
        return pids

    def _kill_pids(self, pids):
        for pid in pids:
            if pid == -1:
                continue
            try:
                os.kill(pid, 9)
            except:
                self.logger.info("PID {0} is unassigned".format(pid))

    def clean_before_serving(self):
        self.logger.debug("Cleaning processes before serving")
        tmp_file = tmp_pid_file(self.model_id)
        dir_name = os.path.dirname(tmp_file)
        pids = []
        for proc_file in os.listdir(dir_name):
            if proc_file[-3:] != "pid":
                continue
            proc_file = os.path.join(dir_name, proc_file)
            self.logger.debug(proc_file)
            pids += self._pids_from_file(proc_file)
            os.remove(proc_file)
        self.logger.debug("Cleaning {0} processes".format(pids))
        self._kill_pids(pids)

    def clean_temp_dir(self):
        self.logger.debug("Cleaning temp dir")
        tmp_folder = tempfile.gettempdir()
        for d in os.listdir(tmp_folder):
            if "ersilia-" in d:
                d = os.path.join(tmp_folder, d)
                self.logger.debug("Flushing temporary directory {0}".format(d))
                try:
                    shutil.rmtree(d)
                except:
                    self.logger.warning(
                        "Could not remove temporary directory {0}".format(d)
                    )

    def clean_docker_containers(self):
        self.logger.debug("Silencing docker containers if necessary")
        dm = DockerManager(config_json=self.config_json)
        if dm.is_inside_docker():
            self.logger.debug("It is inside docker")
            return
        if dm.is_installed():
            self.logger.debug("It is not inside docker")
            dm.stop_containers(self.model_id)

    def serve(self):
        self.clean_before_serving()
        self.clean_temp_dir()
        self.close()
        self.service.serve()
        tmp_file = tmp_pid_file(self.model_id)
        with open(tmp_file, "a+") as f:
            f.write("{0} {1}{2}".format(self.service.pid, self.service.url, os.linesep))

    def close(self):
        tmp_file = tmp_pid_file(self.model_id)
        if os.path.isfile(tmp_file):
            pids = self._pids_from_file(tmp_file)
            self._kill_pids(pids)
            os.remove(tmp_file)
        self.clean_temp_dir()
        self.clean_docker_containers()
        try:
            self.service.close()
        except:  # TODO: capture the error
            pass

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
        _api = Api(
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
