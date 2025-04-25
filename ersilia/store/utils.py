import csv
import gzip
import io
import os
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import click
import requests
from tqdm import tqdm

from ersilia.default import INFERENCE_STORE_API_URL


def log_prefix():
    now = datetime.now().strftime("%H:%M:%S")
    return f"[INFO]:[{now}] "


def echo_intro(click_iface):
    width = 80
    click_iface.echo("╔" + "═" * width + "╗", fg="blue", bold=True)
    click_iface.echo(
        f"║{'Ersilia Precalculation Store Engine':^{width}}║", fg="blue", bold=True
    )
    click_iface.echo(f"║{'Version 1.0':^{width}}║", fg="blue", bold=True)
    click_iface.echo("╚" + "═" * width + "╝", fg="blue", bold=True)


def echo_uploading_inputs(click_iface):
    click_iface.echo(
        f"{log_prefix()}Uploading input CSV to S3...", fg="green", bold=True
    )


def echo_upload_complete(click_iface):
    click_iface.echo(f"{log_prefix()}Upload complete.", fg="green", bold=True)


def echo_submitting_job(click_iface):
    click_iface.echo(
        f"{log_prefix()}Submitting a task for precalculation", fg="blue", bold=True
    )


def echo_job_submitted(click_iface, job_id: str):
    click_iface.echo(
        f"{log_prefix()}Job submitted with a job id: {job_id}", fg="blue", bold=True
    )


def echo_polling_status(click_iface):
    click_iface.echo(
        "Waiting for precalculation to finish. Please be patient!",
        fg="white",
        bg="blue",
    )


def echo_status(click_iface, status: str):
    if status == JobStatus.PENDING:
        click_iface.echo(
            f"{log_prefix()}Submitted job status: {status}", fg="yellow", bold=True
        )
    elif status == JobStatus.FAILED:
        click_iface.echo(
            f"{log_prefix()}Submitted job status: {status}", fg="red", bold=True
        )
    elif status == JobStatus.RUNNING:
        click_iface.echo(
            f"{log_prefix()}Submitted job status: {status}", fg="blue", bold=True
        )
    else:
        click_iface.echo(
            f"{log_prefix()}Submitted job status: {status}", fg="green", bold=True
        )


def echo_job_succeeded(click_iface):
    click_iface.echo(
        f"{log_prefix()}Precalculation successfully fetched", fg="blue", bold=True
    )


def echo_found_shards(click_iface, count: int):
    click_iface.echo(
        f"{log_prefix()}Fetched {count} .gz compressed precalcultaion files.",
        fg="blue",
        bold=True,
    )


def echo_merged_saved(click_iface, output_path: Path):
    click_iface.echo(
        f"{log_prefix()}Merged CSV saved to: {output_path}", fg="blue", bold=True
    )


class ClickInterface:
    """
    Interactive CLI helper for styled console output and progress bars.

    Parameters
    ----------
    colourize : bool, optional
        Whether to enable ANSI styling in output (default is True).
    """

    def __init__(self, colourize: bool = True):
        self.colourize = colourize

    def echo(
        self,
        message: str,
        fg: str = None,
        bg: str = None,
        bold: bool = False,
        blink: bool = False,
    ):
        """
        Print a styled message with optional colors and text effects.

        Parameters
        ----------
        message : str
            The text to display.
        fg : str, optional
            Foreground color name or code.
        bg : str, optional
            Background color name or code.
        bold : bool, optional
            Whether to use bold text.
        blink : bool, optional
            Whether to use blinking text.
        """
        style_args = {}
        if fg:
            style_args["fg"] = fg
        if bg:
            style_args["bg"] = bg
        if bold:
            style_args["bold"] = True
        if blink:
            style_args["blink"] = True
        styled = click.style(message, **style_args) if self.colourize else message
        click.echo(styled)

    def progress_bar(self, total: int, **kwargs):
        """
        Create and return a progress bar instance.

        Parameters
        ----------
        total : int
            Total count for the progress bar.
        **kwargs
            Additional keyword arguments passed to tqdm.

        Returns
        -------
        tqdm.tqdm
            Configured progress bar.
        """
        return tqdm(total=total, **kwargs)


class ApiClient:
    """
    HTTP client wrapper for JSON-based and streaming API requests.

    Parameters
    ----------
    base_url : str, optional
        Base URL for all HTTP requests.
    timeout : int, optional
        Default timeout in seconds for HTTP requests.
    """

    def __init__(self, base_url: str = None, timeout: int = 3600):
        self.base_url = base_url or ""
        self.timeout = timeout

    def get_json(self, endpoint: str, params: dict = None):
        """
        Perform a GET request and return parsed JSON.

        Parameters
        ----------
        endpoint : str
            API endpoint path (appended to base_url).
        params : dict, optional
            Query parameters for the request.

        Returns
        -------
        dict
            Parsed JSON response.
        """
        r = requests.get(endpoint, params=params, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def post_json(self, endpoint: str, json: dict = None):
        """
        Perform a POST request with JSON payload and return parsed JSON.

        Parameters
        ----------
        endpoint : str
            API endpoint path (appended to base_url).
        json : dict, optional
            JSON payload to send.

        Returns
        -------
        dict
            Parsed JSON response.
        """
        r = requests.post(endpoint, json=json, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def head(self, url: str):
        """
        Perform a HEAD request and return response headers.

        Parameters
        ----------
        url : str
            Full URL to send HEAD request to.

        Returns
        -------
        requests.structures.CaseInsensitiveDict
            Response headers.
        """
        r = requests.head(url, timeout=self.timeout)
        r.raise_for_status()
        return r.headers

    def download_stream(self, url: str, chunk_size: int = 64 * 1024):
        """
        Stream binary data from a URL in configurable chunk sizes.

        Parameters
        ----------
        url : str
            URL to download data from.
        chunk_size : int, optional
            Size of each chunk in bytes (default is 65536).

        Yields
        ------
        bytes
            Chunks of data from the response.
        """
        r = requests.get(url, stream=True, timeout=self.timeout)
        r.raise_for_status()
        for chunk in r.iter_content(chunk_size):
            yield chunk


class FileManager:
    """
    Utility class for temporary CSV creation, S3 uploads, and merging result shards.
    """

    @staticmethod
    def create_temp_csv(records, input_adapter) -> str:
        """
        Write records to a temporary CSV file via a GenericInputAdapter.

        Parameters
        ----------
        records : list
            Iterable of input records to write.
        input_adapter : GenericInputAdapter
            Adapter that yields dicts with 'input' key.

        Returns
        -------
        str
            Path of the created temporary CSV file.
        """
        tmp = tempfile.NamedTemporaryFile(
            delete=False, mode="w+", newline="", encoding="utf-8"
        )
        writer = csv.writer(tmp)
        for rec in input_adapter.adapt_one_by_one(records):
            writer.writerow([rec["input"]])
        tmp.flush()
        tmp.close()
        return tmp.name

    @staticmethod
    def upload_to_s3(presigned_info: dict, file_path: str):
        """
        Upload a local file to S3 using provided presigned URL and form fields.

        Parameters
        ----------
        presigned_info : dict
            Contains 'url' and optional 'fields' for multipart upload.
        file_path : str
            Local path of the file to upload.
        """
        with open(file_path, "rb") as fbin:
            files = {"file": (file_path, fbin)}
            up = requests.post(
                presigned_info["url"],
                data=presigned_info.get("fields", {}),
                files=files,
                timeout=3600,
            )
        delete_file_upon_upload(up.status_code, file_path)
        up.raise_for_status()

    @staticmethod
    def merge_shards(
        shards: list,
        header: list,
        output_path: Path,
        downloader,
        click_iface: ClickInterface,
    ):
        """
        Download, decompress, and concatenate multiple CSV shards into a single file.

        Parameters
        ----------
        shards : list
            List of dicts with keys 'url', 'key', and optional 'size'.
        header : list
            CSV header row.
        output_path : Path
            Path to write the merged CSV.
        downloader : ApiClient
            Used for HEAD and streaming GET requests.
        click_iface : ClickInterface
            CLI interface for printing progress and warnings.
        """
        total_size = 0
        for shard in shards:
            size = shard.get("size", 0)
            try:
                headers = downloader.head(shard["url"])
                size = int(headers.get("Content-Length", size))
            except requests.HTTPError:
                pass
            total_size += size
        click_iface.echo(
            f"{log_prefix()}Total download file size: {total_size / 1024:.1f} KB",
            fg="blue",
            bold=True,
        )

        click_iface.echo(
            f"{log_prefix()}Downloading and merging each file content",
            fg="blue",
            bold=True,
        )
        with open(output_path, "w", encoding="utf-8", newline="") as out_f:
            out_f.write(",".join(header) + "\n")
            first = True
            pbar = click_iface.progress_bar(
                total=total_size,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc="⬇ GZIP",
                colour="blue",
                ncols=70,
                bar_format="{desc}: [{bar}] {percentage:3.0f}% {n_fmt}/{total_fmt} • {rate_fmt} • ETA: {remaining}",
            )
            for shard in shards:
                buf = io.BytesIO()
                for chunk in downloader.download_stream(shard["url"]):
                    buf.write(chunk)
                    pbar.update(len(chunk))
                data = buf.getvalue()
                if shard.get("key", "").endswith(".gz") or shard["url"].endswith(".gz"):
                    data = gzip.decompress(data)
                text = data.decode("utf-8", errors="replace")
                for i, line in enumerate(text.splitlines(keepends=True)):
                    if not first and i == 0 and line.startswith(header[0]):
                        continue
                    out_f.write(line)
                first = False
            pbar.close()


class InferenceStoreMessage(object):
    """
    Base class for inference store messages.

    Parameters
    ----------
    model_id : str
        The ID of the model for which the message is being generated.
    """

    def __init__(self, model_id):
        self.model_id = model_id

    def _echo(self, text, **styles):
        return click.echo(click.style(text, **styles))


class JobStatus:
    """
    Class to define job completetion status.
    """

    PENDING = "PENDING"
    RUNNING = "RUNNING"
    FAILED = "FAILED"
    SUCCEEDED = "SUCCEEDED"


class OutputSource:
    """
    Class to define output source options.
    """

    LOCAL_ONLY = "local-cache-only"
    LOCAL = "local-cache"
    CLOUD_ONLY = "cloud-cache-only"
    CLOUD = "cloud-cache"
    ALL = [
        LOCAL_ONLY,
        CLOUD_ONLY,
    ]

    @classmethod
    def is_local(cls, option):
        """
        Check if the option is local.

        Parameters
        ----------
        option : str
            The option to check.

        Returns
        -------
        bool
            True if the option is local, False otherwise.
        """
        return option == cls.LOCAL_ONLY

    @classmethod
    def is_cloud(cls, option):
        """
        Check if the option is cloud.

        Parameters
        ----------
        option : str
            The option to check.

        Returns
        -------
        bool
            True if the option is cloud, False otherwise.
        """
        return option == cls.CLOUD_ONLY


class ModelNotInStore(InferenceStoreMessage):
    """
    Message class for models not found in the inference store.

    Parameters
    ----------
    model_id : str
        The ID of the model that is not found.
    """

    def __init__(self, model_id):
        super().__init__(model_id)
        self.model_id = model_id

    def echo(self):
        """
        Echo the message for model not found in inference store.
        """
        super()._echo(
            "Model {0} could not be found in inference store".format(self.model_id),
            fg="red",
        )
        super()._echo(
            "Please serve the model locally: ersilia serve {0} --output-source {1}".format(
                self.model_id, OutputSource.LOCAL_ONLY
            )
        )
        sys.exit(0)


class PrecalculationsNotInStore(InferenceStoreMessage):
    """
    Message class for precalculations not found in the inference store.

    Parameters
    ----------
    model_id : str
        The ID of the model for which precalculations are not found.
    """

    def __init__(self, model_id):
        super().__init__(model_id)
        self.model_id = model_id

    def echo(self):
        """
        Echo the message for precalculations not found in inference store.
        """
        super()._echo(
            "Precalculations for model {0} could not be found in inference store".format(
                self.model_id
            ),
            fg="red",
        )
        super()._echo(
            "Please serve the model locally: ersilia serve {0} --output-source {1}".format(
                self.model_id, OutputSource.LOCAL_ONLY
            )
        )
        sys.exit(0)


class PrecalculationsInStore(InferenceStoreMessage):
    """
    Message class for precalculations found in the inference store.

    Parameters
    ----------
    model_id : str
        The ID of the model for which precalculations are found.
    output_url : str
        The URL where the precalculations can be downloaded.
    """

    def __init__(self, model_id, output_url):
        super().__init__(model_id)
        self.output_url = output_url

    def echo(self):
        """
        Echo the message for precalculations available for download.

        Parameters
        ----------
        output_url : str
            The URL for downloading the precalculations.
        """
        super()._echo(
            "Precalculations for model {0} are now available for download via this link (expires in 60 minutes): {1}".format(
                self.model_id, self.output_url
            ),
            fg="green",
        )
        sys.exit(0)


def store_has_model(model_id: str) -> bool:
    """
    Check if the model exists in the inference store.

    Parameters
    ----------
    model_id : str
        The ID of the model to check.

    Returns
    -------
    bool
        True if the model exists in the store, False otherwise.
    """
    response = requests.get(
        INFERENCE_STORE_API_URL + "/model", params={"modelid": model_id}, timeout=60
    )
    if response.status_code == 200:
        print(f"Model {model_id} found in inference store")
        return True
    print(f"Model {model_id} not found in inference store")
    return False


def delete_file_upon_upload(response_code: int, file_path: str):
    """
    Delete the file upon successful upload.

    Parameters
    ----------
    response_code : int
        The HTTP response code from the upload request.
    file_path : str
        The path of the file to delete.
    """
    if response_code == 200:
        try:
            os.remove(file_path)
            print(f"File {file_path} deleted successfully.")
        except Exception as e:
            print(
                f"Failed to delete file {file_path}, please delete manually. Error: {e}"
            )
