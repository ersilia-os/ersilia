import json
import os
import subprocess
import threading
import time

import docker
from dockerfile_parse import DockerfileParser

from .. import logger
from ..default import (
    DEFAULT_DOCKER_PLATFORM,
    DEFAULT_UDOCKER_USERNAME,
    DOCKER_INFO_FILE,
    DOCKERHUB_LATEST_TAG,
    DOCKERHUB_ORG,
)
from ..utils.logging import make_temp_dir
from ..utils.system import SystemChecker
from .identifiers.long import LongIdentifier
from .terminal import run_command, run_command_check_output


def set_docker_host():
    try:
        # Get the current Docker context
        context_result = subprocess.run(
            ["docker", "context", "show"], capture_output=True, text=True, check=True
        )
        context_name = context_result.stdout.strip()

        # Get the Docker host for the current context
        result = subprocess.run(
            [
                "docker",
                "context",
                "inspect",
                context_name,
                "--format",
                "{{.Endpoints.docker.Host}}",
            ],
            capture_output=True,
            text=True,
            check=True,
        )

        base_url = result.stdout.strip()

        if base_url:
            os.environ["DOCKER_HOST"] = base_url  # Set the variable for this process
    except:
        return


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


def model_image_version_reader(dir):
    """
    Read the requested model image version from a file.
    """
    if os.path.exists(os.path.join(dir, DOCKER_INFO_FILE)):
        with open(os.path.join(dir, DOCKER_INFO_FILE), "r") as f:
            data = json.load(f)
        if "tag" in data:
            return data["tag"]
    return DOCKERHUB_LATEST_TAG


class SimpleDocker(object):
    """
    A class to manage Docker containers and images.

    Parameters
    ----------
    use_udocker : bool, optional
        Whether to use udocker instead of Docker. Default is None.
    """

    def __init__(self, use_udocker=None):
        self.identifier = LongIdentifier()
        self.logger = logger
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
        """
        Get a dictionary of Docker images.

        Returns
        -------
        dict
            A dictionary of Docker images with image names as keys and image IDs as values.
        """
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
        """
        Get a dictionary of Docker containers.

        Parameters
        ----------
        only_run : bool
            Whether to include only running containers.

        Returns
        -------
        dict
            A dictionary of Docker containers with container names as keys and image names as values.
        """
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
            # cnt_idx = h.find("CONTAINER ID")
            img_idx = h.find("IMAGE")
            cmd_idx = h.find("COMMAND")
            # sts_idx = h.find("STATUS")
            # pts_idx = h.find("PORTS")
            nam_idx = h.find("NAMES")
            for l in f:
                # cnt = l[cnt_idx:img_idx].strip()
                img = l[img_idx:cmd_idx].strip()
                # sts = l[sts_idx:pts_idx].strip()
                nam = l[nam_idx:].strip()
                cnt_dict[nam] = img
        return cnt_dict

    def exists(self, org, img, tag):
        """
        Check if a Docker image exists.

        Parameters
        ----------
        org : str
            The organization name.
        img : str
            The image name.
        tag : str
            The image tag.

        Returns
        -------
        bool
            True if the image exists, False otherwise.
        """
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
        """
        Build a Docker image.

        Parameters
        ----------
        path : str
            The path to the Dockerfile.
        org : str
            The organization name.
        img : str
            The image name.
        tag : str
            The image tag.

        """
        if self._with_udocker:
            raise Exception("Cannot built with udocker")
        path = os.path.abspath(path)
        cwd = os.getcwd()
        os.chdir(path)
        cmd = "docker build -t %s %s" % (self._image_name(org, img, tag), path)
        run_command(cmd)
        os.chdir(cwd)

    def delete(self, org, img, tag):
        """
        Delete a Docker image.

        Parameters
        ----------
        org : str
            The organization name.
        img : str
            The image name.
        tag : str
            The image tag.
        """
        if not self._with_udocker:
            cmd = "docker rmi -f %s" % self._image_name(org, img, tag)
            run_command(cmd)
        else:
            cmd = "sudo -u {0} udocker rmi {1}".format(
                DEFAULT_UDOCKER_USERNAME, self._image_name(org, img, tag)
            )

    def run(self, org, img, tag, name, memory=None):
        """
        Run a Docker container.

        Parameters
        ----------
        org : str
            The organization name.
        img : str
            The image name.
        tag : str
            The image tag.
        name : str
            The container name.
        memory : int, optional
            The memory limit for the container in GB. Default is None.

        Returns
        -------
        str
            The name of the running container.
        """
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
            cmd = "sudo -u {0} udocker run {2} bash".format(DEFAULT_UDOCKER_USERNAME)  # noqa: F524
            run_command(cmd)
        return name

    @staticmethod
    def kill(name):
        """
        Kill a Docker container.

        Parameters
        ----------
        name : str
            The container name.
        """
        cmd = "docker kill {0}".format(name)
        run_command(cmd)

    @staticmethod
    def remove(name):
        """
        Remove a Docker container.

        Parameters
        ----------
        name : str
            The container name.
        """
        cmd = "docker rm -f {0}".format(name)
        run_command(cmd)

    @staticmethod
    def cp_from_container(name, img_path, local_path, org=None, img=None, tag=None):
        """
        Copy files from a Docker container to the local filesystem.

        Parameters
        ----------
        name : str
            The container name.
        img_path : str
            The path inside the container.
        local_path : str
            The local path to copy files to.
        org : str, optional
            The organization name. Default is None.
        img : str, optional
            The image name. Default is None.
        tag : str, optional
            The image tag. Default is None.
        None
        """
        local_path = os.path.abspath(local_path)
        dirname = os.path.dirname(local_path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
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

    async def cp_from_image(self, img_path, local_path, org, img, tag):
        """
        Copy files from a Docker image to the local filesystem.

        Parameters
        ----------
        img_path : str
            The path inside the image.
        local_path : str
            The local path to copy files to.
        org : str
            The organization name.
        img : str
            The image name.
        tag : str
            The image tag.
        """
        name = self.run(org, img, tag, name=None)
        self.cp_from_container(name, img_path, local_path, org=org, img=img, tag=tag)
        self.remove(name)

    @staticmethod
    def exec_container(name, cmd):
        """
        Execute a command in a Docker container.

        Parameters
        ----------
        name : str
            The container name.
        cmd : str
            The command to execute.
        """
        cmd = 'docker exec -i %s bash -c "%s"' % (name, cmd)
        run_command(cmd)

    def exec(self, cmd, org, img, tag, name):
        """
        Execute a command in a Docker container and then kill the container.

        Parameters
        ----------
        cmd : str
            The command to execute.
        org : str
            The organization name.
        img : str
            The image name.
        tag : str
            The image tag.
        name : str
            The container name.
        """
        name = self.run(org, img, tag, name=name)
        self.exec_container(name, cmd)
        self.kill(name)

    def container_peak(self, model_id):
        """
        Get the peak memory usage of a Docker container running an Ersilia model.

        Parameters
        ----------
        model_id : str
            The model identifier.

        Returns
        -------
        float or None
            The peak memory usage in MB, or None if the container is not found or an error occurs.
        """
        try:
            set_docker_host()
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
                    self.logger.debug(
                        f"Could not compute container peak memory for model {model_id}"
                    )
                    return
            else:
                self.logger.debug(f"No container found for model {model_id}")
                return

        except docker.errors.NotFound as e:
            self.logger.debug(f"Container {container.name} not found: {e}")
            return None
        except docker.errors.APIError as e:
            logger.debug(f"Docker API error: {e}")
            return None
        except Exception as e:
            self.logger.debug(f"An error occurred: {e}")
            return None

    def cleanup_ersilia_images(self):
        """
        Remove all Ersilia-related Docker images.
        """
        if self._with_udocker:
            self.logger.warning("Docker cleanup not supported with udocker")
            return

        try:
            images_dict = self.images()

            if not images_dict:
                logger.info("No Docker images found")
                return

            for image_name, image_id in images_dict.items():
                if DOCKERHUB_ORG in image_name:
                    try:
                        logger.info(f"Removing Docker image: {image_name}")
                        self.delete(*self._splitter(image_name))
                    except Exception as e:
                        logger.error(f"Failed to remove Docker image {image_name}: {e}")

        except Exception as e:
            self.logger.error(f"Failed to cleanup Docker images: {e}")


class SimpleDockerfileParser(DockerfileParser):
    """
    A class to parse Dockerfiles.

    Parameters
    ----------
    path : str
        The path to the Dockerfile or the directory containing the Dockerfile.
    """

    def __init__(self, path):
        if os.path.isdir(path):
            path = os.path.join(path, "Dockerfile")
        DockerfileParser.__init__(self, path=path)

    def get_baseimage(self):
        """
        Get the base image from the Dockerfile.

        Returns
        -------
        str
            The base image.
        """
        return self.baseimage

    def get_runs(self):
        """
        Get the RUN commands from the Dockerfile.

        Returns
        -------
        list
            A list of RUN commands.
        """
        structure = self.structure
        runs = []
        for d in structure:
            if d["instruction"] == "RUN":
                val = d["value"]
                for v in val.split("&&"):
                    runs += [v.strip()]
        return runs


class ContainerMetricsSampler:
    """
    A class to sample metrics from Docker containers.

    Parameters
    ----------
    model_id : str
        The model identifier.
    sampling_interval : float, optional
        The interval between samples in seconds. Default is 0.01.
    """

    def __init__(self, model_id, sampling_interval=0.01):
        set_docker_host()
        self.client = docker.from_env()
        self.logger = logger
        self.container = self._get_container_for_model(model_id)
        self.cpu_samples = []
        self.memory_samples = []
        self.tracking_thread = None
        self.tracking = False
        self.sampling_interval = sampling_interval  # Sample every 10ms

    def _get_container_for_model(self, model_id=None):
        """
        This function will get the Docker container running an ersilia model

        Parameters
        ----------
        model_id : str
            The model identifier.

        Returns
        -------
        container
            The Docker container running the model.
        """
        if not model_id:
            raise ValueError("Model ID is required")
        try:
            set_docker_host()
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

        Parameters
        ----------
        stats : dict
            The stats dictionary from the Docker API.

        Returns
        -------
        float
            The CPU percentage utilisation.
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

        Parameters
        ----------
        stats : dict
            The stats dictionary from the Docker API.

        Returns
        -------
        float
            The memory usage percentage.
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
        """
        Start tracking metrics from the Docker container.
        """
        if not self.tracking:
            self.tracking = True
            self.cpu_samples.clear()
            self.memory_samples.clear()
            self.tracking_thread = threading.Thread(target=self._collect_metrics)
            self.tracking_thread.start()

    def stop_tracking(self):
        """
        Stop tracking metrics from the Docker container.
        """
        if self.tracking:
            self.tracking = False
            if self.tracking_thread is not None:
                self.tracking_thread.join()
                self.tracking_thread = None

    def get_average_metrics(self):
        """
        Get the average metrics from the tracked samples.

        Returns
        -------
        dict
            A dictionary containing the average CPU and memory usage, and the peak CPU and memory usage.
        """
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
