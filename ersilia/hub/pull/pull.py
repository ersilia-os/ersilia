import asyncio
import os
import re
import subprocess

import requests

from ... import ErsiliaBase, throw_ersilia_exception
from ...default import DOCKERHUB_LATEST_TAG, DOCKERHUB_ORG
from ...utils.docker import SimpleDocker
from ...utils.echo import echo
from ...utils.exceptions_utils.pull_exceptions import (
    DockerConventionalPullError,
    DockerImageNotAvailableError,
)
from ...utils.logging import make_temp_dir
from ...utils.terminal import run_command, yes_no_input

# Re-downloading an image the user already has is opt-in.
PULL_IMAGE = os.environ.get("PULL_IMAGE", "n")


class PullProgress:
    """
    Progress of a Docker image pull, from the Docker Engine API's pull events.

    The Engine API reports progress as structured events (layer id, status,
    and byte counts), the same across Docker versions, instead of the text
    that ``docker pull`` prints.

    Parameters
    ----------
    expected_bytes : int, optional
        Compressed size of the image, if known beforehand. Used as the total
        until every layer has reported its own size.
    layer_sizes : dict, optional
        Compressed size of each layer, keyed by its short digest (the 12
        characters Docker uses as layer id), if known beforehand (e.g. from
        Docker Hub). Docker only reports a layer's size once it starts
        downloading, so this keeps the total exact from the start.
    """

    DONE = ("Pull complete", "Already exists")

    def __init__(self, expected_bytes=None, layer_sizes=None):
        self.expected_bytes = expected_bytes or 0
        self.layer_sizes = layer_sizes or {}
        self.layers = {}

    def update(self, event):
        """
        Record one pull event.

        Parameters
        ----------
        event : dict
            A decoded event, e.g. ``{"id": "55d2dadd4bbc", "status":
            "Downloading", "progressDetail": {"current": 1, "total": 2}}``.
        """
        layer_id = event.get("id")
        status = event.get("status") or ""
        # Events without a layer (e.g. "Pulling from ...", "Digest: ...") and
        # the first event, whose id is the tag, carry no layer progress.
        if not layer_id or status.startswith("Pulling from"):
            return
        detail = event.get("progressDetail") or {}
        layer = self.layers.setdefault(
            layer_id,
            {
                "size": self.layer_sizes.get(layer_id, 0),
                "downloaded": 0,
                "extracted": 0,
                "status": "",
            },
        )
        layer["status"] = status
        total = detail.get("total") or 0
        current = detail.get("current") or 0
        if total:
            layer["size"] = max(layer["size"], total)
        if status.startswith("Downloading"):
            layer["downloaded"] = current
        elif status in ("Verifying Checksum", "Download complete"):
            layer["downloaded"] = layer["size"]
        elif status.startswith("Extracting"):
            layer["downloaded"] = layer["size"]
            layer["extracted"] = current
        elif status == "Pull complete":
            layer["downloaded"] = layer["extracted"] = layer["size"]

    @property
    def layers_done(self):
        """Number of layers pulled or already present."""
        return sum(1 for layer in self.layers.values() if layer["status"] in self.DONE)

    @property
    def layers_total(self):
        """Number of layers seen so far."""
        return len(self.layers)

    @property
    def total_bytes(self):
        """Bytes to download (compressed), as far as known."""

        def exists(layer_id):
            # Layers already present locally are not downloaded.
            layer = self.layers.get(layer_id)
            return layer is not None and layer["status"] == "Already exists"

        if self.layer_sizes:
            # Every listed layer counts from the start, until Docker reports
            # it as already present; so the total can only shrink, early on.
            listed = sum(
                size
                for layer_id, size in self.layer_sizes.items()
                if not exists(layer_id)
            )
            unlisted = sum(
                layer["size"]
                for layer_id, layer in self.layers.items()
                if layer_id not in self.layer_sizes and not exists(layer_id)
            )
            return listed + unlisted
        present = [
            layer for layer_id, layer in self.layers.items() if not exists(layer_id)
        ]
        known = sum(layer["size"] for layer in present)
        if len(present) < len(self.layers):
            return known
        return max(known, self.expected_bytes)

    @property
    def downloaded_bytes(self):
        """Bytes downloaded so far."""
        return sum(layer["downloaded"] for layer in self.layers.values())

    @property
    def extracted_bytes(self):
        """Bytes extracted so far."""
        return sum(layer["extracted"] for layer in self.layers.values())

    @property
    def extracting(self):
        """Whether everything is downloaded and layers are being extracted."""
        # Docker starts extracting a layer while others are still downloading;
        # the text only switches to "extracting" once downloading is done.
        return self.downloaded_bytes >= self.total_bytes and any(
            layer["status"].startswith("Extracting") for layer in self.layers.values()
        )

    def bar(self):
        """
        Completed and total units for a progress bar.

        Returns
        -------
        tuple of (float, float or None)
            Downloading and extracting count as two halves when byte counts
            are known; otherwise finished layers out of all layers. The total
            is None until anything is known.
        """
        total = self.total_bytes
        if total:
            return self.downloaded_bytes + self.extracted_bytes, 2 * total
        if self.layers_total:
            return self.layers_done, self.layers_total
        return 0, None

    def describe(self):
        """
        Text shown next to the bar, e.g. "95/182 MB  3/14 layers".

        Returns
        -------
        str
            Megabytes (decimal, as Docker reports them) and layers.
        """
        parts = []
        total = self.total_bytes
        if total:
            done = self.extracted_bytes if self.extracting else self.downloaded_bytes
            prefix = "extracting " if self.extracting else ""
            parts.append(f"{prefix}{done / 1e6:.0f}/{total / 1e6:.0f} MB")
        if self.layers_total:
            parts.append(f"{self.layers_done}/{self.layers_total} layers")
        return "  ".join(parts)


def pull_with_progress(
    repository,
    tag,
    callback,
    platform=None,
    expected_bytes=None,
    layer_sizes=None,
    stop=None,
):
    """
    Pull a Docker image through the Docker Engine API, reporting progress.

    Parameters
    ----------
    repository : str
        Image repository, e.g. "ersiliaos/eos3b5e".
    tag : str
        Image tag.
    callback : callable
        Called with a ``PullProgress`` after every event.
    platform : str, optional
        Platform to pull, e.g. "linux/amd64".
    expected_bytes : int, optional
        Compressed image size, if known.
    layer_sizes : dict, optional
        Compressed size of each layer by short digest, if known.
    stop : threading.Event, optional
        When set, the pull stops (e.g. after Ctrl+C).

    Raises
    ------
    RuntimeError
        If the pull fails.
    """
    import docker

    from ...utils.docker import set_docker_host

    set_docker_host()
    client = docker.from_env(timeout=600)
    progress = PullProgress(expected_bytes=expected_bytes, layer_sizes=layer_sizes)
    try:
        events = client.api.pull(
            repository, tag=tag, stream=True, decode=True, platform=platform
        )
        for event in events:
            if stop is not None and stop.is_set():
                return
            if event.get("error"):
                raise RuntimeError(event["error"])
            progress.update(event)
            callback(progress)
    except docker.errors.APIError as e:
        raise RuntimeError(str(e)) from e
    except (requests.exceptions.RequestException, docker.errors.DockerException) as e:
        raise RuntimeError("connection to Docker lost: {0}".format(e)) from e


def pull_error(model_id, tag, error):
    """
    Turn a pull failure into an error that says why and what to do.

    Parameters
    ----------
    model_id : str
        The model identifier.
    tag : str
        The image tag that was pulled.
    error : Exception
        The failure reported by Docker.

    Returns
    -------
    ImagePullError
        The error to raise.
    """
    from ...utils.exceptions_utils.cli_exceptions import ImagePullError

    text = str(error).lower()
    if "no space left on device" in text:
        return ImagePullError(
            model_id,
            "Docker ran out of disk space",
            "Free space with 'docker system prune', or raise the disk limit in Docker Desktop > Settings > Resources.",
        )
    if "toomanyrequests" in text or "rate limit" in text:
        return ImagePullError(
            model_id,
            "the Docker Hub download limit was reached",
            "Log in with 'docker login', or try again later.",
        )
    if "manifest unknown" in text or "not found" in text:
        return ImagePullError(
            model_id,
            "version '{0}' does not exist".format(tag),
            "See the available versions at https://hub.docker.com/r/{0}/{1}/tags".format(
                DOCKERHUB_ORG, model_id
            ),
        )
    if "connection to docker lost" in text:
        return ImagePullError(
            model_id,
            "Docker stopped responding during the download",
            "Start Docker and run the fetch again.",
        )
    return ImagePullError(model_id, str(error), "")


class ModelPuller(ErsiliaBase):
    """
    ModelPuller is responsible for pulling models from DockerHub.

    Parameters
    ----------
    model_id : str
        The ID of the model to be pulled.
    overwrite : bool, optional
        Whether to overwrite existing files.
    config_json : dict, optional
        Configuration settings for the puller.

    Examples
    --------
    .. code-block:: python

        puller = ModelPuller(
            model_id="eosxxxx", config_json=config
        )
        await puller.async_pull()
    """

    def __init__(
        self,
        model_id: str,
        overwrite: bool = None,
        config_json: dict = None,
        docker_tag: str = None,
    ):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.simple_docker = SimpleDocker()
        self.model_id = model_id
        self.docker_tag = docker_tag or DOCKERHUB_LATEST_TAG
        self.image_name = "{0}/{1}:{2}".format(
            DOCKERHUB_ORG, self.model_id, self.docker_tag
        )
        self.overwrite = overwrite

    def is_available_locally(self) -> bool:
        """
        Check if the Docker image is available locally.

        Returns
        -------
        bool
            True if the image is available locally, False otherwise.
        """
        is_available = self.simple_docker.exists(
            DOCKERHUB_ORG, self.model_id, self.docker_tag
        )
        if is_available:
            self.logger.debug("Image {0} is available locally".format(self.image_name))
            return True
        else:
            self.logger.debug(
                "Image {0} is not available locally".format(self.image_name)
            )
            return False

    def is_available_in_dockerhub(self) -> bool:
        """
        Check if the Docker image is available in DockerHub.

        Returns
        -------
        bool
            True if the image is available in DockerHub, False otherwise.
        """
        from ...utils.exceptions_utils.cli_exceptions import ImagePullError

        url = "https://hub.docker.com/v2/repositories/{0}/{1}/tags/{2}".format(
            DOCKERHUB_ORG, self.model_id, self.docker_tag
        )
        try:
            response = requests.get(url, timeout=15)
        except requests.exceptions.RequestException as e:
            raise ImagePullError(
                self.model_id,
                "Docker Hub could not be reached",
                "Check your internet connection and try again.",
            ) from e
        if response.status_code == 429:
            raise pull_error(self.model_id, self.docker_tag, "toomanyrequests")
        if response.status_code == 404 and self.docker_tag != DOCKERHUB_LATEST_TAG:
            raise pull_error(self.model_id, self.docker_tag, "manifest unknown")
        if response.status_code == 200:
            self.logger.debug(
                "The docker image {0} exists in DockerHub".format(self.image_name)
            )
            return True
        else:
            self.logger.debug(
                "The docker image {0} does not exist in DockerHub".format(
                    self.image_name
                )
            )
            return False

    def _delete(self):
        self.logger.debug(
            "Deleting locally available image {0}".format(self.image_name)
        )
        self.simple_docker.delete(
            org=DOCKERHUB_ORG, img=self.model_id, tag=self.docker_tag
        )

    def _get_size_of_local_docker_image_in_mb(self) -> float:
        try:
            image_name = "{0}/{1}:{2}".format(
                DOCKERHUB_ORG, self.model_id, self.docker_tag
            )
            result = subprocess.check_output(
                ["docker", "image", "inspect", image_name, "--format", "{{.Size}}"]
            )
            size_in_mb = int(result.strip()) / (1024 * 1024)
            return size_in_mb
        except subprocess.CalledProcessError:
            self.logger.warning("Image not found locally")
            return None

    @staticmethod
    def _architecture(platform=None):
        # The architecture Docker pulls: the given platform, or this machine's.
        if platform:
            return platform.split("/")[-1]
        import platform as _platform

        machine = _platform.machine().lower()
        return "arm64" if machine in ("arm64", "aarch64") else "amd64"

    def _get_remote_image_info(self, platform=None):
        # Compressed size and per-layer sizes of the image Docker will pull,
        # from Docker Hub. (None, {}) if unknown.
        url = "https://hub.docker.com/v2/repositories/{0}/{1}/tags/{2}/images".format(
            DOCKERHUB_ORG, self.model_id, self.docker_tag
        )
        try:
            response = requests.get(url, timeout=10)
            if response.status_code != 200:
                return None, {}
            arch = self._architecture(platform)
            for image in response.json():
                if image.get("architecture") != arch:
                    continue
                layers = {
                    layer["digest"].split(":")[-1][:12]: layer.get("size") or 0
                    for layer in image.get("layers") or []
                    if layer.get("digest")
                }
                return image.get("size") or None, layers
        except Exception:
            pass
        return None, {}

    def _get_remote_image_size_mb(self, platform=None) -> float:
        size, _ = self._get_remote_image_info(platform)
        return size / 1e6 if size else None

    @throw_ersilia_exception()
    async def async_pull(self):
        """
        Asynchronously pull the Docker image.
        """
        if self.is_available_locally():
            if self.overwrite is None:
                do_pull = yes_no_input(
                    "The image of model {0} is already available locally. Download it again?".format(
                        self.model_id
                    ),
                    default_answer=PULL_IMAGE,
                )
            elif self.overwrite:
                do_pull = True
            else:
                do_pull = False
            if not do_pull:
                self.logger.info("Skipping pulling the image")
                return
            # The local image is not deleted first: 'docker pull' only
            # replaces it once the new download has succeeded.
        else:
            self.logger.debug("Docker image of the model is not available locally")
        if self.is_available_in_dockerhub():
            self.logger.debug(
                "Pulling image {0} from DockerHub...".format(self.image_name)
            )

            verbose = getattr(self.logger, "verbosity", 0) == 1

            remote_bytes, layer_sizes = (
                self._get_remote_image_info() if not verbose else (None, {})
            )
            remote_size = remote_bytes / 1e6 if remote_bytes else None
            if not verbose:
                size_text = (
                    f" (~{remote_size:.0f} MB compressed)" if remote_size else ""
                )
                echo(f"Downloading the Docker image{size_text}.")

            pull_command = (
                f"docker pull {DOCKERHUB_ORG}/{self.model_id}:{self.docker_tag}"
            )
            force_pull_command = f"docker pull {DOCKERHUB_ORG}/{self.model_id}:{self.docker_tag} --platform linux/amd64"

            if verbose:

                async def _run_pull(cmd):
                    proc = await asyncio.create_subprocess_shell(
                        cmd,
                        stdout=asyncio.subprocess.PIPE,
                        stderr=asyncio.subprocess.PIPE,
                    )

                    async def log_stream(stream, log_method):
                        async for line in stream:
                            log_method(line.decode().strip())

                    await asyncio.gather(
                        log_stream(proc.stdout, self.logger.info),
                        log_stream(proc.stderr, self.logger.error),
                    )
                    await proc.wait()
                    if proc.returncode != 0:
                        raise subprocess.CalledProcessError(proc.returncode, cmd)

                try:
                    await _run_pull(pull_command)
                except subprocess.CalledProcessError:
                    self.logger.warning("Conventional pull failed, trying linux/amd64")
                    await _run_pull(force_pull_command)
            else:
                from rich.progress import (
                    BarColumn,
                    Progress,
                    TextColumn,
                    TimeElapsedColumn,
                )

                repository = f"{DOCKERHUB_ORG}/{self.model_id}"

                with Progress(
                    # Indented to line up under the text of the line above.
                    TextColumn("    "),
                    BarColumn(),
                    TextColumn("{task.fields[detail]}"),
                    TimeElapsedColumn(),
                ) as progress:
                    task = progress.add_task("", total=None, detail="")

                    def on_progress(pull):
                        completed, total = pull.bar()
                        progress.update(
                            task,
                            completed=completed,
                            total=total,
                            detail=pull.describe(),
                        )

                    import threading

                    loop = asyncio.get_running_loop()
                    stop = threading.Event()

                    async def in_thread(platform):
                        # A daemon thread, stopped on Ctrl+C, so an interrupt
                        # does not wait for the whole download to finish.
                        done = loop.create_future()

                        def work():
                            try:
                                result = pull(platform)
                            except BaseException as e:
                                loop.call_soon_threadsafe(done.set_exception, e)
                            else:
                                loop.call_soon_threadsafe(done.set_result, result)

                        threading.Thread(target=work, daemon=True).start()
                        try:
                            return await done
                        except BaseException:
                            stop.set()
                            raise

                    def pull(platform):
                        if platform:
                            expected, sizes = self._get_remote_image_info(platform)
                        else:
                            expected, sizes = remote_bytes, layer_sizes
                        return pull_with_progress(
                            repository,
                            self.docker_tag,
                            on_progress,
                            platform=platform,
                            expected_bytes=expected,
                            layer_sizes=sizes,
                            stop=stop,
                        )

                    try:
                        await in_thread(None)
                    except (KeyboardInterrupt, asyncio.CancelledError) as e:
                        e.ersilia_note = (
                            "Download cancelled. Run the fetch again to resume."
                        )
                        raise
                    except RuntimeError as e:
                        if "no matching manifest" not in str(e).lower():
                            raise pull_error(self.model_id, self.docker_tag, e) from e
                        # No image for this machine's architecture.
                        self.logger.warning(f"Pull failed ({e}), trying linux/amd64")
                        echo(
                            "No image for this machine's architecture; using the Intel (amd64) image, which runs slower.",
                            fg="yellow",
                        )
                        try:
                            await in_thread("linux/amd64")
                        except RuntimeError as e:
                            raise pull_error(self.model_id, self.docker_tag, e) from e
                    # Show the bar full once Docker reports the pull as done.
                    total = progress.tasks[0].total
                    if total:
                        progress.update(task, completed=total)

                echo("Docker image downloaded.", fg="green")

            size = self._get_size_of_local_docker_image_in_mb()
            if size:
                self.logger.debug("Size of image {0} MB".format(size))
            else:
                self.logger.warning("Could not obtain size of image")
            self.simple_docker.label_with_current_user(
                DOCKERHUB_ORG, self.model_id, self.docker_tag
            )
            return size
        else:
            self.logger.info("Image {0} is not available".format(self.image_name))
            raise DockerImageNotAvailableError(model=self.model_id)

    @throw_ersilia_exception()
    def pull(self):
        """
        This method pulls the Docker image non-asynchronously.
        """
        if self.is_available_locally():
            if self.overwrite is None:
                do_pull = yes_no_input(
                    "The image of model {0} is already available locally. Download it again?".format(
                        self.model_id
                    ),
                    default_answer=PULL_IMAGE,
                )
            elif self.overwrite:
                do_pull = True
            else:
                do_pull = False
            if not do_pull:
                self.logger.info("Skipping pulling the image")
                return
            # The local image is not deleted first: 'docker pull' only
            # replaces it once the new download has succeeded.
        else:
            self.logger.debug("Docker image of the model is not available locally")
        if self.is_available_in_dockerhub():
            self.logger.debug(
                "Pulling image {0} from DockerHub...".format(self.image_name)
            )
            try:
                self.logger.debug(
                    "Trying to pull image {0}/{1}".format(DOCKERHUB_ORG, self.model_id)
                )
                tmp_file = os.path.join(
                    make_temp_dir(prefix="ersilia-"), "docker_pull.log"
                )
                self.logger.debug("Keeping logs of pull in {0}".format(tmp_file))
                run_command(
                    "docker pull {0}/{1}:{2} > {3} 2>&1".format(
                        DOCKERHUB_ORG, self.model_id, self.docker_tag, tmp_file
                    )
                )
                with open(tmp_file, "r") as f:
                    pull_log = f.read()
                    self.logger.debug(pull_log)
                if re.search(r"no match.*manifest", pull_log):
                    self.logger.warning(
                        "No matching manifest for image {0}".format(self.model_id)
                    )
                    raise DockerConventionalPullError(model=self.model_id)
                self.logger.debug("Image pulled succesfully!")
            except DockerConventionalPullError:
                self.logger.warning(
                    "Conventional pull did not work, Ersilia is now forcing linux/amd64 architecture"
                )
                run_command(
                    "docker pull {0}/{1}:{2} --platform linux/amd64".format(
                        DOCKERHUB_ORG, self.model_id, self.docker_tag
                    )
                )
            size = self._get_size_of_local_docker_image_in_mb()
            if size:
                self.logger.debug("Size of image {0} MB".format(size))
                # path = os.path.join(self._model_path(self.model_id), MODEL_SIZE_FILE)
                # with open(path, "w") as f:
                #     json.dump({"size": size, "units": "MB"}, f, indent=4)
                # self.logger.debug("Size written to {}".format(path))
            else:
                self.logger.warning("Could not obtain size of image")
            self.simple_docker.label_with_current_user(
                DOCKERHUB_ORG, self.model_id, self.docker_tag
            )
            return size
        else:
            self.logger.info("Image {0} is not available".format(self.image_name))
            raise DockerImageNotAvailableError(model=self.model_id)
