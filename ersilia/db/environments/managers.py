import os
import re
import shutil
import sys

from ...core.base import ErsiliaBase
from ...default import DOCKERHUB_LATEST_TAG, DOCKERHUB_ORG
from ...setup.requirements.docker import DockerRequirement
from ...utils.docker import SimpleDocker, model_image_version_reader, resolve_platform
from ...utils.identifiers.short import ShortIdentifier
from ...utils.logging import make_temp_dir
from ...utils.paths import Paths
from ...utils.ports import find_free_port
from ...utils.session import get_session_dir
from ...utils.system import is_inside_docker
from ...utils.terminal import run_command
from .localdb import EnvironmentDb

BENTOML_DOCKERPORT = 5000
INTERNAL_DOCKERPORT = 80


class DockerManager(ErsiliaBase):
    """
    Manages Docker operations for Ersilia models.

    It provides methods to build, run, and manage Docker images and containers
    associated with Ersilia models.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for initializing the Docker manager.
    preferred_port : int, optional
        Preferred port number for running Docker containers.
    with_bentoml : bool, optional
        Flag indicating whether to use BentoML for model serving.

    Examples
    --------
    >>> docker_manager = DockerManager(
    ...     config_json=config, preferred_port=8080
    ... )
    >>> docker_manager.build(
    ...     model_id="eosxxxx",
    ...     docker_user="user",
    ...     docker_pwd="pass",
    ... )
    >>> docker_manager.run(model_id="eosxxxx", workers=2)
    """

    def __init__(self, config_json=None, preferred_port=None, with_bentoml=False):
        ErsiliaBase.__init__(self, config_json=config_json)
        self._eos_regex = Paths()._eos_regex()
        self._org_regex = re.compile(DOCKERHUB_ORG)
        self.inside_docker = is_inside_docker()
        self.docker = SimpleDocker()
        self.db = EnvironmentDb()
        self.db.table = "docker"
        self.preferred_port = preferred_port
        self.with_bentoml = with_bentoml

    def is_inside_docker(self):
        """
        Checks if the current process is running inside a Docker container.

        Returns
        -------
        bool
            True if running inside Docker, False otherwise.
        """
        return self.inside_docker

    def is_installed(self):
        """
        Checks if Docker is installed on the system.

        Returns
        -------
        bool
            True if Docker is installed, False otherwise.
        """
        return DockerRequirement().is_installed()

    def is_logged_in(self):
        """
        Checks if the user is logged in to Docker.

        Returns
        -------
        bool
            True if the user is logged in to Docker, False otherwise.
        """
        return DockerRequirement().is_logged_in()

    def is_active(self):
        """
        Checks if Docker is active and running.

        Returns
        -------
        bool
            True if Docker is active, False otherwise.
        """
        return DockerRequirement().is_active()

    def image_exists(self, model_id):
        """
        Checks if a Docker image for the specified model exists.

        Parameters
        ----------
        model_id : str
            Identifier of the model.

        Returns
        -------
        bool
            True if the Docker image exists, False otherwise.
        """
        if self.images_of_model(model_id, only_latest=True):
            return True
        else:
            return False

    def container_exists(self, container_name):
        """
        Checks if a Docker container with the specified name exists.

        Parameters
        ----------
        container_name : str
            Name of the Docker container.

        Returns
        -------
        bool
            True if the container exists, False otherwise.
        """
        names = set(self.containers(only_run=False).keys())
        if container_name in names:
            return True
        else:
            return False

    def images(self):
        """
        Retrieves a dictionary of Docker images related to Ersilia models.

        Returns
        -------
        dict
            Dictionary of Docker images with image names as keys.
        """
        img_dict = {}
        for k, v in self.docker.images().items():
            if self._eos_regex.search(k) and self._org_regex.search(k):
                img_dict[k] = v
        return img_dict

    def containers(self, only_run):
        """
        Retrieves a dictionary of Docker containers related to Ersilia models.

        Parameters
        ----------
        only_run : bool
            If True, only return running containers.

        Returns
        -------
        dict
            Dictionary of Docker containers with container names as keys.
        """
        images = self.images()
        cnt_dict = {}
        for k, v in self.docker.containers(only_run=only_run).items():
            if v in images:
                cnt_dict[k] = v
        return cnt_dict

    def images_of_model(self, model_id, only_latest=True):
        """
        Retrieves Docker images for a specific model.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        only_latest : bool, optional
            If True, only return images with the latest tag.

        Returns
        -------
        dict
            Dictionary of Docker images for the model.
        """
        images = self.images()
        img_dict = {}
        for k, v in images.items():
            _, img, tag = self.docker._splitter(k)
            if only_latest:
                if tag != "latest":
                    continue
            if img == model_id:
                img_dict[k] = v
        return img_dict

    def containers_of_model(self, model_id, only_run, only_latest=True):
        """
        Retrieves Docker containers for a specific model.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        only_run : bool
            If True, only return running containers.
        only_latest : bool, optional
            If True, only consider images with the latest tag.

        Returns
        -------
        dict
            Dictionary of Docker containers for the model.
        """
        valid_images = set()
        for k, v in self.images_of_model(
            model_id=model_id, only_latest=only_latest
        ).items():
            valid_images.update([k])
        cnt_dict = {}
        for k, v in self.containers(only_run=only_run).items():
            if v in valid_images:
                cnt_dict[k] = v
        return cnt_dict

    def build_with_bentoml(self, model_id, use_cache=True):  # Ignore for versioning
        """
        Builds a Docker image for the model using BentoML.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        use_cache : bool, optional
            If True, use Docker's cache when building the image.
        """
        bundle_path = self._get_bundle_location(model_id)
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "build.sh")
        cmdlines = ["cd {0}".format(bundle_path)]
        if use_cache:
            cache_str = ""
        else:
            cache_str = "--disable-local-cache"
        cmdlines += [
            "docker build {0} -t {1}/{2}:{3} .".format(
                cache_str, DOCKERHUB_ORG, model_id, DOCKERHUB_LATEST_TAG
            )
        ]
        self.logger.debug(cmdlines)
        with open(tmp_file, "w") as f:
            for cmd in cmdlines:
                f.write(cmd + os.linesep)
        cmd = "bash {0}".format(tmp_file)
        run_command(cmd)
        self.db.insert(
            model_id=model_id,
            env="{0}/{1}:{2}".format(DOCKERHUB_ORG, model_id, DOCKERHUB_LATEST_TAG),
        )

    @property
    def _model_deploy_dockerfiles_url(self):
        return "https://raw.githubusercontent.com/ersilia-os/ersilia/master/dockerfiles/model-deploy"

    def _build_ersilia_base(self):
        self.logger.debug("Creating docker image of ersilia base")
        path = make_temp_dir(prefix="ersilia-")
        base_folder = os.path.join(path, "base")
        os.mkdir(base_folder)
        base_files = [
            "Dockerfile",
            "docker-entrypoint.sh",
            "nginx.conf",
            "download-miniconda.sh",
        ]
        for f in base_files:
            cmd = "cd {0}; wget {2}/base/{1}".format(
                base_folder, f, self._model_deploy_dockerfiles_url
            )
            run_command(cmd)
        cmd = "cd {0}; docker build -t {1}/base:{2} .".format(
            base_folder, DOCKERHUB_ORG, DOCKERHUB_LATEST_TAG
        )
        run_command(cmd)

    def build_with_ersilia(
        self, model_id, docker_user, docker_pwd
    ):  # Ignore for versioning
        """
        Builds a Docker image for the model using Ersilia's base image.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        docker_user : str
            Docker Hub username.
        docker_pwd : str
            Docker Hub password.
        """
        self.logger.debug("Creating docker image with ersilia incorporated")
        if self.image_exists("base"):
            pass
        else:
            self._build_ersilia_base()
        path = make_temp_dir(prefix="ersilia-model")
        model_folder = os.path.join(path, model_id)
        os.mkdir(model_folder)
        cmd = "cd {0}; wget {1}/model/Dockerfile".format(
            model_folder, self._model_deploy_dockerfiles_url
        )
        run_command(cmd)
        file_path = os.path.join(model_folder, "Dockerfile")
        with open(file_path, "r") as f:
            text = f.read()
        text = text.replace("eos_identifier", model_id)
        with open(file_path, "w") as f:
            f.write(text)
        # login to docker so as to be able to push images to dockerhub
        cmd = "docker login --password {0} --username {1}".format(
            docker_pwd, docker_user
        )
        run_command(cmd)

        try:
            # Attempt to build for linux/amd64 and linux/arm64
            print("building for both linux/amd64,linux/arm64")
            cmd = "cd {0}; docker buildx build --builder=container --platform linux/amd64,linux/arm64 -t {1}/{2}:{3} --push .".format(
                model_folder, DOCKERHUB_ORG, model_id, DOCKERHUB_LATEST_TAG
            )
            run_command(cmd)
            print("done building for linux/amd64,linux/arm64")
            sys.exit()  # This will terminate the program immediately only if build for linux/amd64,linux/arm64 complete successfully

        except:
            # Build failed for multi-platforms, now try building only for linux/amd64
            self.logger.warning(
                "Build failed for multi-platform, trying linux/amd64 only"
            )
            cmd = "cd {0}; docker build -t {1}/{2}:{3} .".format(
                model_folder, DOCKERHUB_ORG, model_id, DOCKERHUB_LATEST_TAG
            )
            run_command(cmd)

    def build(self, model_id, docker_user, docker_pwd, use_cache=True):
        """
        Builds a Docker image for the model.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        docker_user : str
            Docker Hub username.
        docker_pwd : str
            Docker Hub password.
        use_cache : bool, optional
            If True, use Docker's cache when building the image.
        """
        if self.with_bentoml:
            self.build_with_bentoml(model_id=model_id, use_cache=use_cache)
        else:
            self.build_with_ersilia(
                model_id=model_id, docker_user=docker_user, docker_pwd=docker_pwd
            )

    def remove(self, model_id):
        """
        Removes the Docker image for the specified model.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        """
        bundle_path = self._get_bundle_location(model_id)
        docker_tag = model_image_version_reader(bundle_path)
        self.docker.delete(org=DOCKERHUB_ORG, img=model_id, tag=docker_tag)
        self.db.delete(
            model_id=model_id,
            env="{0}/{1}:{2}".format(DOCKERHUB_ORG, model_id, docker_tag),
        )

    def run(self, model_id, workers=1, enable_microbatch=True, memory=None):
        """
        Runs a Docker container for the specified model.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        workers : int, optional
            Number of worker processes.
        enable_microbatch : bool, optional
            If True, enable micro-batching.
        memory : int, optional
            Memory limit for the container in gigabytes.

        Returns
        -------
        dict
            Dictionary containing container name and port.
        """
        self.logger.debug("Running docker manager")
        imgs = self.images_of_model(model_id, only_latest=True)
        self.logger.debug("Available images: {0}".format(imgs))
        img = [k for k in imgs.keys()][0]
        if enable_microbatch:
            mb_string = "--enable-microbatch"
        else:
            mb_string = ""
        port = find_free_port(preferred_port=self.preferred_port)
        name = None
        si = ShortIdentifier()
        while not name:
            name_ = "{0}_{1}".format(model_id, si.encode())
            if not self.container_exists(name_):
                name = name_
        if self.with_bentoml:
            dockerport = BENTOML_DOCKERPORT
        else:
            dockerport = INTERNAL_DOCKERPORT
        if memory is None:
            cmd = "docker run --platform {6} --name {0} -d -p {1}:{2} {3} --workers={4} {5}".format(
                name, port, dockerport, img, workers, mb_string, resolve_platform()
            )
        else:
            cmd = 'docker run --memory="{7}" --platform {6} --name {0} -d -p {1}:{2} {3} --workers={4} {5}'.format(
                name,
                port,
                dockerport,
                img,
                workers,
                mb_string,
                resolve_platform(),
                str(memory) + "g",
            )
        self.logger.debug(cmd)
        run_command(cmd)
        return {"container_name": name, "port": port}

    def start(self, container_name):
        """
        Starts a stopped Docker container.

        Parameters
        ----------
        container_name : str
            Name of the Docker container to start.
        """
        cmd = "docker start {0}".format(container_name)
        run_command(cmd)

    def stop(self, container_name):
        """
        Stops a running Docker container.

        Parameters
        ----------
        container_name : str
            Name of the Docker container to stop.
        """
        cmd = "docker stop {0}".format(container_name)
        run_command(cmd)

    def _delete_container(self, container_name):
        cmd = "docker rm --force {0}".format(container_name)
        run_command(cmd)

    def delete_containers(self, model_id):
        """
        Deletes all Docker containers associated with a model.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        """
        for k, _ in self.containers_of_model(
            model_id, only_run=False, only_latest=False
        ).items():
            self._delete_container(k)

    def delete_image(self, model_id):
        """
        Deletes a Docker image associated with a model.
        """
        self.remove(model_id)

    def remove_stopped_containers(self):
        """
        Removes all stopped Docker containers.
        """
        if self.is_inside_docker():
            return
        if not self.is_installed():
            return
        cmd = "docker container prune -f"
        self.logger.debug("Removing stopped containers")
        run_command(cmd)

    def _stop_containers_with_model_id(self, model_id):
        if self.is_inside_docker():
            return
        if not self.is_installed():
            return
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "docker-ps.txt")
        cmd = "docker ps > {0}".format(tmp_file)
        self.logger.debug("Running {0}".format(cmd))
        run_command(cmd)
        cids = []
        with open(tmp_file, "r") as f:
            h = next(f)
            col_idx = len(h.split("IMAGE")[0])
            for l in f:
                img_str = l[col_idx:].split(" ")[0]
                if model_id in img_str:
                    cid = l.split(" ")[0]
                    cids += [cid]
        for cid in cids:
            cmd = "docker container kill {0}".format(cid)
            run_command(cmd)
        shutil.rmtree(tmp_folder)

    def _stop_containers_with_entrypoint_sh(self):
        if self.is_inside_docker:
            return
        if not self.is_installed():
            return
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "docker-ps.txt")
        cmd = "docker ps > {0}".format(tmp_file)
        self.logger.debug("Running {0}".format(cmd))
        run_command(cmd)
        cids = []
        with open(tmp_file, "r") as f:
            h = next(f)
            col_idx = len(h.split("COMMAND")[0])
            for l in f:
                cmd_str = l[col_idx:].split(" ")[0]
                if "entrypoint.sh bash" in cmd_str:
                    cid = l.split(" ")[0]
                    cids += [cid]
        for cid in cids:
            cmd = "docker container kill {0}".format(cid)
            run_command(cmd)
        shutil.rmtree(tmp_folder)

    def stop_containers(self, model_id):
        """
        Stops all running Docker containers associated with a model.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        """
        if self.is_inside_docker:
            return
        if not self.is_installed():
            return
        self._stop_containers_with_model_id(model_id)
        self._stop_containers_with_entrypoint_sh()
        self.remove_stopped_containers()
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "docker-ps.txt")
        cmd = "docker ps > {0}".format(tmp_file)
        self.logger.debug("Running {0}".format(cmd))
        run_command(cmd)
        cids = []
        with open(tmp_file, "r") as f:
            next(f)
            for l in f:
                mid = l.rstrip().split(" ")[-1].split("_")[0]
                if mid == model_id:
                    cid = l.split(" ")[0]
                    cids += [cid]
        for cid in cids:
            cmd = "docker container kill {0}".format(cid)
            run_command(cmd)
        shutil.rmtree(tmp_folder)
        self.remove_stopped_containers()

    def prune(self):
        """
        Prunes unused Docker objects to free up space.
        """
        cmd = "docker system prune -f"
        run_command(cmd)

    def delete_image(self, img):  # noqa: D102, F811
        fn = os.path.join(get_session_dir(), "rm_image_output.txt")
        cmd = "docker image rm {0} --force 2> {1}".format(img, fn)
        run_command(cmd)
        with open(fn, "r") as f:
            text = f.read()
            patt = "image is being used by running container "
            if patt in text:
                container_id = text.split(patt)[1].rstrip()
                self.logger.debug(
                    "A running container was found {0}. Removing it before the image".format(
                        container_id
                    )
                )
                cmd = "docker stop {0}".format(container_id)
                run_command(cmd)
                cmd = "docker rm {0}".format(container_id)
                run_command(cmd)
                self.delete_image(img)

    def delete_images(self, model_id, purge_unnamed=True):
        """
        Deletes Docker images associated with a model.

        Parameters
        ----------
        model_id : str
            Identifier of the model.
        purge_unnamed : bool, optional
            If True, also remove unnamed images.
        """
        if self.is_inside_docker():
            return
        if not self.is_installed():
            return
        self.stop_containers(model_id)
        self.prune()
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "docker-images.txt")
        cmd = "docker images > {0}".format(tmp_file)
        self.logger.debug("Running {0}".format(cmd))
        run_command(cmd)
        unnamed_images = []
        named_images = []
        with open(tmp_file, "r") as f:
            h = next(f)
            img_idx = len(h.split("IMAGE ID")[0])
            for l in f:
                img = l[img_idx:].split(" ")[0]
                name = l.split(" ")[0]
                if model_id in name:
                    named_images += [img]
                if "<none>" in name:
                    unnamed_images += [img]
        images = named_images
        if purge_unnamed:
            images += unnamed_images
        for img in images:
            self.logger.debug("Removing docker image {0}".format(img))
            self.delete_image(img)


class CondaManager(object):  # noqa: D101
    def __init__(self):
        pass

    def environments(self):  # noqa: D102
        pass
