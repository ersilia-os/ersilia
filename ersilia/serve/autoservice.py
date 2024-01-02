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


class AutoService(ErsiliaBase):
    def __init__(
        self,
        model_id,
        service_class=None,
        config_json=None,
        preferred_port=None,
        url=None,
    ):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.logger.debug("Setting AutoService for {0}".format(model_id))
        self.config_json = config_json
        self.model_id = model_id
        self._meta = None
        self._preferred_port = preferred_port
        self._url = url
        if service_class is None:
            self.logger.debug("No service class provided, deciding automatically")
            service_class_file = os.path.join(
                self._get_bundle_location(model_id), SERVICE_CLASS_FILE
            )
            if os.path.exists(service_class_file):
                self.logger.debug(
                    "Service class file exists in folder {0}".format(service_class_file)
                )
                with open(service_class_file, "r") as f:
                    s = f.read()
                if not s:
                    s = None
            else:
                s = None
            if s is not None:
                self.logger.debug("Service class: {0}".format(s))
                if s == "system":
                    self.service = SystemBundleService(
                        model_id, config_json=config_json, preferred_port=preferred_port
                    )
                elif s == "venv":
                    self.service = VenvEnvironmentService(
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
            else:
                self.logger.debug(
                    "No service class file exists in {0}".format(service_class_file)
                )
                with open(service_class_file, "w") as f:
                    if SystemBundleService(
                        model_id, config_json=config_json, preferred_port=preferred_port
                    ).is_available():
                        self.service = SystemBundleService(
                            model_id,
                            config_json=config_json,
                            preferred_port=preferred_port,
                        )
                        self.logger.debug("Service class: system")
                        f.write("system")
                        self._service_class = "system"
                    elif VenvEnvironmentService(
                        model_id, config_json=config_json, preferred_port=preferred_port
                    ).is_available():
                        self.service = VenvEnvironmentService(
                            model_id,
                            config_json=config_json,
                            preferred_port=preferred_port,
                        )
                        f.write("venv")
                        self.logger.debug("Service class: venv")
                        self._service_class = "venv"
                    elif CondaEnvironmentService(
                        model_id, config_json=config_json, preferred_port=preferred_port
                    ).is_available():
                        self.service = CondaEnvironmentService(
                            model_id,
                            config_json=config_json,
                            preferred_port=preferred_port,
                        )
                        f.write("conda")
                        self.logger.debug("Service class: conda")
                        self._service_class = "conda"
                    elif DockerImageService(
                        model_id, config_json=config_json, preferred_port=preferred_port
                    ).is_available():
                        self.service = DockerImageService(
                            model_id,
                            config_json=config_json,
                            preferred_port=preferred_port,
                        )
                        f.write("docker")
                        self.logger.debug("Service class: docker")
                        self._service_class = "docker"
                    elif PulledDockerImageService(
                        model_id, config_json=config_json, preferred_port=preferred_port
                    ).is_available():
                        self.service = PulledDockerImageService(
                            model_id,
                            config_json=config_json,
                            preferred_port=preferred_port,
                        )
                        f.write("pulled_docker")
                        self.logger.debug("Service class: pulled_docker")
                        self._service_class = "pulled_docker"
                    elif HostedService(
                        model_id, config_json=config_json, url=url
                    ).is_available():
                        self.service = HostedService(
                            model_id, config_json=config_json, url=url
                        )
                        f.write("hosted")
                        self.logger.debug("Service class: hosted")
                        self._service_class = "hosted"
                    else:
                        self.logger.debug("Service class: dummy")
                        self.service = DummyService(
                            model_id,
                            config_json=config_json,
                            preferred_port=preferred_port,
                        )
        else:
            self.logger.info("Service class provided")
            service_class = self._service_class_loader(service_class)
            if service_class(
                model_id,
                config_json=config_json,
                preferred_port=preferred_port,
                url=url,
            ).is_available():
                self.service = service_class(
                    model_id,
                    config_json=config_json,
                    preferred_port=preferred_port,
                    url=url,
                )
            else:
                self.service = None
        self._set_apis()

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

    def _service_class_loader(self, service_class):
        if type(service_class) is SystemBundleService:
            self._service_class = "system"
            return service_class
        elif type(service_class) is VenvEnvironmentService:
            self._service_class = "venv"
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
            elif service_class == "venv":
                return VenvEnvironmentService
            elif service_class == "conda":
                return CondaEnvironmentService
            elif service_class == "docker":
                return DockerImageService
            elif service_class == "pulled_docker":
                return PulledDockerImageService
            elif service_class == "hosted":
                return HostedService
            raise Exception()

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
            if api_name not in self._meta:
                self._meta = {api_name: _api.meta()}
            yield result
