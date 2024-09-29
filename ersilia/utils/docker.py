import os
import docker
import subprocess
import threading
import time
import json
from dockerfile_parse import DockerfileParser

from .identifiers.long import LongIdentifier
from .terminal import run_command, run_command_check_output

from .. import logger
from ..default import (DEFAULT_DOCKER_PLATFORM, DEFAULT_UDOCKER_USERNAME,
                       DOCKERHUB_ORG, DOCKERHUB_LATEST_TAG,
                       PACK_METHOD_BENTOML, PACK_METHOD_FASTAPI)
from ..utils.system import SystemChecker
from ..utils.logging import make_temp_dir

def resolve_pack_method_docker(model_id):
    client = docker.from_env()
    model_image = client.images.get(
        f"{DOCKERHUB_ORG}/{model_id}:{DOCKERHUB_LATEST_TAG}"
    )
    image_history = model_image.history()
    for hist in image_history:
        # Very hacky, but works bec we don't have nginx in ersilia-pack images
        if "nginx" in hist["CreatedBy"]: 
            return PACK_METHOD_BENTOML
    return PACK_METHOD_FASTAPI


def resolve_platform():
    if SystemChecker().is_arm64():
        return "linux/arm64"
    else:
        return DEFAULT_DOCKER_PLATFORM


def is_docker_installed():
    try:
        subprocess.check_output(["docker", "--version"])
        return True
    except subprocess.CalledProcessError:
        return False
    except FileNotFoundError:
        return False


def is_udocker_installed():
    cmd = ["sudo", "-u", DEFAULT_UDOCKER_USERNAME, "udocker", "--help"]
    try:
        out = subprocess.check_output(cmd)
        if "Syntax" in out:
            return True
        return False
    except:
        return False


class SimpleDocker(object):
    def __init__(self, use_udocker=None):
        self.identifier = LongIdentifier()
        if use_udocker is None:
            self._with_udocker = self._use_udocker()
        else:
            self._with_udocker = use_udocker

    def _use_udocker(self):
        if is_docker_installed():
            return False
        if not is_udocker_installed():
            return False
        return True

    @staticmethod
    def _splitter(name):
        name = name.split("/")
        if len(name) == 2:
            org = name[0]
            img = name[1]
        else:
            org = None
            img = name[0]
        img_ = img.split(":")
        if len(img_) == 2:
            img = img_[0]
            tag = img_[1]
        else:
            img = img_[0]
        return org, img, tag

    @staticmethod
    def _image_name(org, img, tag):
        return "%s/%s:%s" % (org, img, tag)

    def images(self):
        tmp_dir = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_dir, "images.txt")
        if not self._with_udocker:
            cmd = "docker images > {0}".format(tmp_file)
            run_command(cmd)
            img_dict = {}
            with open(tmp_file, "r") as f:
                h = next(f)
                rep_idx = h.find("REPOSITORY")
                tag_idx = h.find("TAG")
                img_idx = h.find("IMAGE ID")
                crt_idx = h.find("CREATED")
                for l in f:
                    rep = l[rep_idx:tag_idx].strip()
                    tag = l[tag_idx:img_idx].strip()
                    img = l[img_idx:crt_idx].strip()
                    img_dict["{0}:{1}".format(rep, tag)] = img
            return img_dict
        else:
            # TODO
            cmd = "sudo -u {0} udocker images > {1}".format(
                DEFAULT_UDOCKER_USERNAME, tmp_file
            )
            run_command(cmd)
            with open(tmp_file, "r") as f:
                h = next(f)
                for l in f:
                    img = l.split(" ")[0]
                    img_dict[img] = img
            return img_dict

    def containers(self, only_run):
        tmp_dir = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_dir, "containers.txt")
        if not only_run:
            all_str = "-a"
        else:
            all_str = ""
        cmd = "docker ps {0} > {1}".format(all_str, tmp_file)
        run_command(cmd)
        cnt_dict = {}
        with open(tmp_file, "r") as f:
            h = next(f)
            cnt_idx = h.find("CONTAINER ID")
            img_idx = h.find("IMAGE")
            cmd_idx = h.find("COMMAND")
            sts_idx = h.find("STATUS")
            pts_idx = h.find("PORTS")
            nam_idx = h.find("NAMES")
            for l in f:
                cnt = l[cnt_idx:img_idx].strip()
                img = l[img_idx:cmd_idx].strip()
                sts = l[sts_idx:pts_idx].strip()
                nam = l[nam_idx:].strip()
                cnt_dict[nam] = img
        return cnt_dict

    def exists(self, org, img, tag):
        bash_script = """
        #!/bin/bash
        image_and_tag="$1"
        image_and_tag_array=(${image_and_tag//:/ })
        if [[ "$(docker images ${image_and_tag_array[0]} | grep ${image_and_tag_array[1]} 2> /dev/null)" != "" ]]; then
          echo "True"
        else
          echo "False"
        fi
        """
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_script = os.path.join(tmp_folder, "exists.sh")
        with open(tmp_script, "w") as f:
            f.write(bash_script)
        cmd = ["bash", tmp_script, self._image_name(org, img, tag)]
        res = run_command_check_output(cmd)
        res = res.strip()
        if res == b"False":
            return False
        if res == b"True":
            return True
        return None

    def build(self, path, org, img, tag):
        if self._with_udocker:
            raise Exception("Cannot built with udocker")
        path = os.path.abspath(path)
        cwd = os.getcwd()
        os.chdir(path)
        cmd = "docker build -t %s %s" % (self._image_name(org, img, tag), path)
        run_command(cmd)
        os.chdir(cwd)

    def delete(self, org, img, tag):
        if not self._with_udocker:
            cmd = "docker rmi -f %s" % self._image_name(org, img, tag)
            run_command(cmd)
        else:
            cmd = "sudo -u {0} udocker rmi {1}".format(
                DEFAULT_UDOCKER_USERNAME, self._image_name(org, img, tag)
            )

    def run(self, org, img, tag, name, memory=None):
        if name is None:
            name = self.identifier.encode()
        if not self._with_udocker:
            if memory is None:
                cmd = "docker run -it -d --platform {0} --name {1} {2} bash".format(
                    resolve_platform(), name, self._image_name(org, img, tag)
                )
            else:
                cmd = 'docker run -it -d --memory="{3}" --platform {0} --name {1} {2} bash'.format(
                    resolve_platform(),
                    name,
                    self._image_name(org, img, tag),
                    str(memory) + "g",
                )
            run_command(cmd)
        else:
            # TODO
            cmd = "sudo -u {0} udocker run {2} bash".format(
                DEFAULT_UDOCKER_USERNAME, self._image_name(org, img, tag)
            )
            run_command(cmd)
        return name

    @staticmethod
    def kill(name):
        cmd = "docker kill {0}".format(name)
        run_command(cmd)

    @staticmethod
    def remove(name):
        cmd = "docker rm -f {0}".format(name)
        run_command(cmd)

    @staticmethod
    def cp_from_container(name, img_path, local_path, org=None, img=None, tag=None):
        local_path = os.path.abspath(local_path)
        tmp_file = os.path.join(make_temp_dir(prefix="ersilia-"), "tmp.txt")
        cmd = "docker cp %s:%s %s &> %s" % (name, img_path, local_path, tmp_file)
        run_command(cmd)
        with open(tmp_file, "r") as f:
            output = f.read()
        if "No such container" in output:
            img_name = "{0}/{1}:{2}".format(org, img, tag)
            cmd = "docker run --platform linux/amd64 -d --name={0} {1}".format(
                name, img_name
            )
            run_command(cmd)
            cmd = "docker cp %s:%s %s" % (name, img_path, local_path)
            run_command(cmd)

    def cp_from_image(self, img_path, local_path, org, img, tag):
        name = self.run(org, img, tag, name=None)
        self.cp_from_container(name, img_path, local_path, org=org, img=img, tag=tag)
        self.remove(name)

    @staticmethod
    def exec_container(name, cmd):
        cmd = 'docker exec -i %s bash -c "%s"' % (name, cmd)
        run_command(cmd)

    def exec(self, cmd, org, img, tag, name):
        name = self.run(org, img, tag, name=name)
        self.exec_container(name, cmd)
        self.kill(name)

    def container_peak(self, model_id):
        """
        This function will get the peak memory of the Docker container running Ersilia Models.
        """
        try:
            client = docker.from_env()
            for ctr in client.containers.list():
                if model_id in ctr.name:
                    container = ctr
                    break
            peak_memory = None
            if container:
                stats = container.stats(stream=False)
                if "memory_stats" in stats and "max_usage" in stats["memory_stats"]:
                    peak_memory = stats["memory_stats"]["max_usage"] / (1024 * 1024)
                    return peak_memory
                if peak_memory is None:
                    cgroup_path = None
                    possible_paths = [
                        f"/sys/fs/cgroup/system.slice/docker-{container.id}.scope/memory.peak",
                        f"/sys/fs/cgroup/docker/{container.id}/memory.peak",
                    ]
                    for path in possible_paths:
                        if os.path.exists(path):
                            cgroup_path = path
                            break

                    if cgroup_path is None:
                        print(
                            f"Could not find cgroup file for container '{container.name}'"
                        )
                    try:
                        with open(cgroup_path, "r") as file:
                            peak_memory = int(file.read().strip()) / (1024 * 1024)
                    except FileNotFoundError:
                        print(f"cgroup file {cgroup_path} not found")
                    except Exception as e:
                        print(
                            f"An error occurred while reading cgroup file {cgroup_path}: {e}"
                        )
                if peak_memory is not None:
                    return peak_memory
                else:
                    logger.debug(
                        f"Could not compute container peak memory for model {model_id}"
                    )
                    return
            else:
                logger.debug(f"No container found for model {model_id}")
                return

        except docker.errors.NotFound as e:
            print(f"Container {container.name} not found: {e}")
            return None
        except docker.errors.APIError as e:
            print(f"Docker API error: {e}")
            return None
        except Exception as e:
            print(f"An error occurred: {e}")
            return None


class SimpleDockerfileParser(DockerfileParser):
    def __init__(self, path):
        if os.path.isdir(path):
            path = os.path.join(path, "Dockerfile")
        DockerfileParser.__init__(self, path=path)

    def get_baseimage(self):
        return self.baseimage

    def get_runs(self):
        structure = self.structure
        runs = []
        for d in structure:
            if d["instruction"] == "RUN":
                val = d["value"]
                for v in val.split("&&"):
                    runs += [v.strip()]
        return runs


class ContainerMetricsSampler:
    #TODO Perhaps this can be moved to a base Sampling class
    def __init__(self, model_id, sampling_interval=0.01):
        self.client = docker.from_env()
        self.logger = logger
        self.container = self._get_container_for_model(model_id)
        self.cpu_samples = []
        self.memory_samples = []
        self.tracking_thread = None
        self.tracking = False
        self.sampling_interval = sampling_interval # Sample every 10ms

    def _get_container_for_model(self, model_id=None):
        """
        This function will get the Docker container running an ersilia model
        """
        if not model_id:
            raise ValueError("Model ID is required")
        try:
            client = docker.from_env()
            containers = client.containers.list()

            if not containers:
                self.logger.debug("No containers found")
            for container in containers:
                if model_id in container.name:
                    return container

        except docker.errors.APIError as e:
            self.logger.debug(f"Error: Docker API error: {e}")
        except KeyError as e:
            self.logger.debug(f"KeyError: {e} in stats for container.")
        except Exception as e:
            self.logger.debug(f"An unexpected error occurred: {e}")

    def _container_cpu(self, stats: dict):
        """
        This function will get the CPU percentage utilisation of the Docker container running Ersilia Models.
        """

        try:  # Ref: https://stackoverflow.com/a/77924494/1887515
            cpu_usage = (
                stats["cpu_stats"]["cpu_usage"]["total_usage"]
                - stats["precpu_stats"]["cpu_usage"]["total_usage"]
            )
            cpu_system = (
                stats["cpu_stats"]["system_cpu_usage"]
                - stats["precpu_stats"]["system_cpu_usage"]
            )
            num_cpus = stats["cpu_stats"]["online_cpus"]
            cpu_perc = round((float(cpu_usage) / float(cpu_system)) * num_cpus * 100.0)
            self.logger.debug(f"CPU Percentage: {cpu_perc}")
            return cpu_perc
        except KeyError:
            self.logger.debug(
                "KeyError: 'cpu_stats' or 'precpu_stats' not found in stats for container."
            )
            return

    def _container_memory(self, stats: dict):
        """
        This function will get the total memory usage of the Docker container running an ersilia model.
        """
        try:
            mem_usage_mb = stats["memory_stats"]["usage"] / (1024 * 1024)
            mem_available_mb = stats["memory_stats"]["limit"] / (1024 * 1024)
            mem_perc = round((float(mem_usage_mb) / float(mem_available_mb)) * 100)
            return mem_perc
        except KeyError:
            self.logger.debug(
                "KeyError: 'memory_stats' not found in stats for container."
            )
            return

    def _collect_metrics(self):
        if self.container:
            while self.tracking:
                stat = self.container.stats(stream=False)
                cpu_usage = self._container_cpu(stat)
                memory_usage = self._container_memory(stat)
                self.cpu_samples.append(cpu_usage)
                self.memory_samples.append(memory_usage)
                time.sleep(self.sampling_interval)  
        else:
            self.logger.debug("No container found for model")

    def start_tracking(self):
        if not self.tracking:
            self.tracking = True
            self.cpu_samples.clear()
            self.memory_samples.clear()
            self.tracking_thread = threading.Thread(target=self._collect_metrics)
            self.tracking_thread.start()

    def stop_tracking(self):
        if self.tracking:
            self.tracking = False
            if self.tracking_thread is not None:
                self.tracking_thread.join()
                self.tracking_thread = None

    def get_average_metrics(self):
        metrics = {
            "container_cpu_perc": 0,
            "container_memory_perc": 0,
            "peak_memory_perc": 0,
            "peak_cpu_perc": 0,
        }
        if not self.cpu_samples or not self.memory_samples:
            return metrics

        metrics["container_cpu_perc"] = sum(self.cpu_samples) / len(self.cpu_samples)
        metrics["container_memory_perc"] = sum(self.memory_samples) / len(
            self.memory_samples
        )
        metrics["peak_memory"] = max(self.memory_samples)
        metrics["peak_cpu_perc"] = max(self.cpu_samples)
        self.logger.debug(f"Average CPU Percentage: {metrics['container_cpu_perc']}")
        self.logger.debug(
            f"Average Memory Percentage: {metrics['container_memory_perc']}"
        )
        return metrics
