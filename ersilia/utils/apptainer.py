import os
import subprocess

import requests

from ..default import SIF_BUCKET_URL, SIF_DIR
from .echo import echo
from .exceptions_utils.cli_exceptions import (
    ApptainerDownloadError,
    ApptainerImageNotFoundError,
)

# Versions are probed in order (v1, v2, ...); stop after this many, in case the
# bucket ever answers 200 for everything.
MAX_VERSIONS = 50
CHUNK_BYTES = 1 << 20

# Where ersilia-pack puts model bundles inside the image. The image's own
# entrypoint points elsewhere, so the bundle is looked up instead.
BUNDLE_LOOKUP = "ls -d ${ERSILIA_PATH:-/opt/ersilia}/bundles/*/ 2>/dev/null | head -1"


def sif_name(model_id, version):
    """
    File name of a model's image, the same in the bucket and on disk.

    Parameters
    ----------
    model_id : str
        The model identifier, e.g. "eos4e40".
    version : str
        The image version, e.g. "v1".

    Returns
    -------
    str
        E.g. "eos4e40_v1.sif".
    """
    return "{0}_{1}.sif".format(model_id, version)


def sif_path(model_id, version):
    """
    Where a model's image is kept on this machine.

    Parameters
    ----------
    model_id : str
        The model identifier.
    version : str
        The image version.

    Returns
    -------
    str
        A path under ``SIF_DIR`` (``~/eos/sifs`` unless ``ERSILIA_SIF_DIR`` is set).
    """
    return os.path.join(SIF_DIR, sif_name(model_id, version))


def sif_url(model_id, version):
    """
    URL of a model's image in the SIF bucket.

    Parameters
    ----------
    model_id : str
        The model identifier.
    version : str
        The image version.

    Returns
    -------
    str
        The URL, built from ``SIF_BUCKET_URL``.
    """
    return "{0}/{1}".format(SIF_BUCKET_URL, sif_name(model_id, version))


def remote_size(model_id, version):
    """
    Size of a model's image in the bucket.

    The bucket cannot be listed, so a missing image answers 403, not 404.

    Parameters
    ----------
    model_id : str
        The model identifier.
    version : str
        The image version.

    Returns
    -------
    int or None
        Size in bytes, or None if there is no such image.

    Raises
    ------
    ApptainerDownloadError
        If the bucket cannot be reached.
    """
    try:
        response = requests.head(sif_url(model_id, version), timeout=30)
    except requests.exceptions.RequestException as e:
        raise ApptainerDownloadError(
            model_id, "the image store could not be reached"
        ) from e
    if response.status_code != 200:
        return None
    try:
        return int(response.headers["Content-Length"])
    except (KeyError, ValueError):
        return None


def latest_version(model_id):
    """
    The newest version of a model's image.

    Versions are "v1", "v2", ... and the bucket cannot be listed, so they are
    tried in order until one is missing.

    Parameters
    ----------
    model_id : str
        The model identifier.

    Returns
    -------
    str
        The newest version, e.g. "v2".

    Raises
    ------
    ApptainerImageNotFoundError
        If the model has no image at all.
    """
    latest = None
    for n in range(1, MAX_VERSIONS + 1):
        version = "v{0}".format(n)
        if remote_size(model_id, version) is None:
            break
        latest = version
    if latest is None:
        raise ApptainerImageNotFoundError(model_id)
    return latest


def download(model_id, version, verbose=False):
    """
    Download a model's image, showing progress like a Docker pull.

    Nothing is downloaded if a complete copy is already on disk. The file is
    written under a temporary name and renamed at the end, so an interrupted
    download never looks complete.

    Parameters
    ----------
    model_id : str
        The model identifier.
    version : str
        The image version.
    verbose : bool, optional
        Log progress instead of drawing a progress bar.

    Returns
    -------
    str
        Path of the image on disk.

    Raises
    ------
    ApptainerImageNotFoundError
        If there is no such image.
    ApptainerDownloadError
        If the download fails or is incomplete.
    """
    expected = remote_size(model_id, version)
    if expected is None:
        raise ApptainerImageNotFoundError(model_id, version)
    path = sif_path(model_id, version)
    if os.path.isfile(path) and os.path.getsize(path) == expected:
        echo("The Apptainer image is already downloaded.")
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    partial = path + ".part"
    echo("Downloading the Apptainer image (~{0:.0f} MB).".format(expected / 1e6))
    try:
        with requests.get(sif_url(model_id, version), stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(partial, "wb") as f:
                _write_with_progress(r, f, expected, verbose)
        if os.path.getsize(partial) != expected:
            raise ApptainerDownloadError(model_id, "the download is incomplete")
        os.replace(partial, path)
    except requests.exceptions.RequestException as e:
        raise ApptainerDownloadError(model_id, str(e)) from e
    finally:
        if os.path.exists(partial):
            os.remove(partial)
    return path


def _write_with_progress(response, handle, expected, verbose):
    chunks = response.iter_content(chunk_size=CHUNK_BYTES)
    if verbose:
        for chunk in chunks:
            handle.write(chunk)
        return
    from rich.progress import BarColumn, Progress, TextColumn, TimeElapsedColumn

    # The same columns as the Docker pull (hub/pull/pull.py).
    with Progress(
        TextColumn("    "),
        BarColumn(),
        TextColumn("{task.fields[detail]}"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("", total=expected, detail="")
        done = 0
        for chunk in chunks:
            handle.write(chunk)
            done += len(chunk)
            progress.update(
                task,
                completed=done,
                detail="{0:.0f}/{1:.0f} MB".format(done / 1e6, expected / 1e6),
            )


class SimpleApptainer(object):
    """
    Runs commands inside an Apptainer image.

    Parameters
    ----------
    binary : str
        The Apptainer command, "apptainer" or "singularity".
    use_unshare : bool, optional
        Prefix commands with ``unshare -r``. Some environments (e.g. Colab,
        which runs as root without a full user namespace setup) need it.
    """

    def __init__(self, binary, use_unshare=False):
        self.binary = binary
        self.use_unshare = use_unshare

    def exec_args(self, sif):
        """
        The command prefix that runs something inside an image.

        Parameters
        ----------
        sif : str
            Path of the image.

        Returns
        -------
        list of str
            E.g. ``["apptainer", "exec", "/path/eos4e40_v1.sif"]``.
        """
        prefix = ["unshare", "-r"] if self.use_unshare else []
        return prefix + [self.binary, "exec", sif]

    def run(self, sif, command, timeout=300):
        """
        Run a shell command inside an image and return what it printed.

        Parameters
        ----------
        sif : str
            Path of the image.
        command : str
            A shell command, run with ``sh -c``.
        timeout : int, optional
            Seconds to wait.

        Returns
        -------
        subprocess.CompletedProcess
            The finished command.
        """
        return subprocess.run(
            self.exec_args(sif) + ["sh", "-c", command],
            capture_output=True,
            text=True,
            timeout=timeout,
        )

    def works(self, sif):
        """
        Checks that commands can run inside the image with these settings.

        Parameters
        ----------
        sif : str
            Path of the image.

        Returns
        -------
        bool
            True if a trivial command succeeds.
        """
        try:
            return self.run(sif, "true", timeout=120).returncode == 0
        except (OSError, subprocess.SubprocessError):
            return False

    def find_bundle(self, sif):
        """
        Find the model bundle inside the image.

        Parameters
        ----------
        sif : str
            Path of the image.

        Returns
        -------
        str or None
            The bundle folder passed to ``ersilia_model_serve --bundle_path``,
            e.g. "/opt/ersilia/bundles/eos4e40", or None if there is none.
        """
        found = self.run(sif, BUNDLE_LOOKUP).stdout.strip().rstrip("/")
        return found or None

    def find_bundle_file(self, sif, bundle, relative_path):
        """
        Find a file in the bundle's version folder inside the image.

        ersilia-pack keeps the model's files one level below the bundle, in a
        folder named after its version (e.g. ``bundles/eos4e40/20250101/``).

        Parameters
        ----------
        sif : str
            Path of the image.
        bundle : str
            The bundle folder, from ``find_bundle``.
        relative_path : str
            Path within the version folder, e.g. "information.json".

        Returns
        -------
        str or None
            Absolute path inside the image, or None if the file does not exist.
        """
        command = "ls -d {0}/*/{1} 2>/dev/null | head -1".format(bundle, relative_path)
        found = self.run(sif, command).stdout.strip()
        return found or None

    def read_file(self, sif, path):
        """
        Read a file inside the image.

        Parameters
        ----------
        sif : str
            Path of the image.
        path : str
            Absolute path of the file inside the image.

        Returns
        -------
        bytes or None
            The contents, or None if the file does not exist.
        """
        result = subprocess.run(
            self.exec_args(sif) + ["cat", path], capture_output=True, timeout=120
        )
        if result.returncode != 0:
            return None
        return result.stdout
