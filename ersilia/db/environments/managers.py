import os
import tempfile
import re
import collections

from ...core.base import ErsiliaBase
from ...default import DOCKERHUB_ORG, DOCKERHUB_LATEST_TAG
from ...utils.paths import Paths
from ...utils.terminal import run_command
from ...utils.docker import SimpleDocker
from ...utils.identifiers.short import ShortIdentifier
from ...utils.ports import find_free_port
from .localdb import EnvironmentDb


BENTOML_DOCKERPORT = 5000


class DockerManager(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self._eos_regex = Paths()._eos_regex()
        self._org_regex = re.compile(DOCKERHUB_ORG)
        self.docker = SimpleDocker()
        self.db = EnvironmentDb()
        self.db.table = "docker"

    def image_exists(self, model_id):
        if self._images_of_model(model_id, only_latest=True):
            return True
        else:
            return False

    def container_exists(self, container_name):
        names = set(self.containers(only_run=False).keys())
        if container_name in names:
            return True
        else:
            return False

    def images(self):
        img_dict = {}
        for k,v in self.docker.images().items():
            if self._eos_regex.search(k) and self._org_regex.search(k):
                img_dict[k] = v
        return img_dict

    def containers(self, only_run):
        images = self.images()
        cnt_dict = {}
        for k,v in self.docker.containers(only_run=only_run).items():
            if v in images:
                cnt_dict[k] = v
        return cnt_dict

    def images_of_model(self, model_id, only_latest=True):
        images = self.images()
        regex = re.compile("{0}/{1}".format(DOCKERHUB_ORG, model_id))
        img_dict = {}
        for k,v in images.items():
            if regex.search(k):
                if only_latest and k.split(":")[-1] != DOCKERHUB_LATEST_TAG:
                    continue
                img_dict[k] = v
        return img_dict

    def containers_of_model(self, model_id, only_run, only_latest=True):
        valid_images = set()
        for k,v in self.images_of_model(model_id=model_id, only_latest=only_latest).items():
            valid_images.update([k])
        cnt_dict = {}
        for k,v in self.containers(only_run=only_run).items():
            if v in valid_images:
                cnt_dict[k] = v
        return cnt_dict

    def build(self, model_id, use_cache=True):
        bundle_path = self._get_bundle_location(model_id)
        tmp_folder = tempfile.mkdtemp()
        tmp_file = os.path.join(tmp_folder, "build.sh")
        cmdlines = ["cd {0}".format(bundle_path)]
        if use_cache:
            cache_str = ""
        else:
            cache_str = "--no-cache"
        cmdlines += ["docker build {0} -t {1}/{2}:{3} .".format(cache_str, DOCKERHUB_ORG, model_id, DOCKERHUB_LATEST_TAG)]
        with open(tmp_file, "w") as f:
            for cmd in cmdlines:
                f.write(cmd+os.linesep)
        cmd = "bash {0}".format(tmp_file)
        run_command(cmd, quiet=False)
        self.db.insert(model_id=model_id, env="{0}/{1}:{2}".format(DOCKERHUB_ORG, model_id, DOCKERHUB_LATEST_TAG))

    def remove(self, model_id):
        self.docker.delete(org=DOCKERHUB_ORG, img=model_id, tag=DOCKERHUB_LATEST_TAG)
        self.db.delete(model_id=model_id, env="{0}/{1}:{2}".format(DOCKERHUB_ORG, model_id, DOCKERHUB_LATEST_TAG))

    def run(self, model_id, workers=1, enable_microbatch=True):
        imgs = self.images_of_model(model_id, only_latest=True)
        img = [k for k in imgs.keys()][0]
        if enable_microbatch:
            mb_string = "--enable-microbatch"
        else:
            mb_string = ""
        port = find_free_port()
        name = None
        si = ShortIdentifier()
        while not name:
            name_ = "{0}_{1}".format(model_id, si.encode())
            if not self.container_exists(name_):
                name = name_
        cmd = "docker run --name {0} -d -p {1}:{2} {3} --workers={4} {5}".format(name, port, BENTOML_DOCKERPORT, img, workers, mb_string)
        run_command(cmd, quiet=True)
        return {"container_name": name, "port": port}

    def start(self, container_name):
        cmd = "docker start {0}".format(container_name)
        run_command(cmd, quiet=True)

    def stop(self, container_name):
        cmd = "docker stop {0}".format(container_name)
        run_command(cmd, quiet=True)

    def _delete_container(self, container_name):
        cmd = "docker rm --force {0}".format(container_name)
        run_command(cmd, quiet=True)

    def delete_containers(self, model_id):
        for k,v in self.containers_of_model(model_id, only_run=False, only_latest=False).items():
            self._delete_container(k)

    def delete_image(self, model_id):
        self.remove(model_id)



class CondaManager(object):

    def __init__(self):
        pass

    def environments(self):
        pass
