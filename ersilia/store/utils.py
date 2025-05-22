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


# ruff: noqa: W291
title = r"""
    ____                            __              __        __   _                  _____  __                   
   / __ \ _____ ___   _____ ____ _ / /_____ __  __ / /____ _ / /_ (_)____   ____     / ___/ / /_ ____   _____ ___ 
  / /_/ // ___// _ \ / ___// __ `// // ___// / / // // __ `// __// // __ \ / __ \    \__ \ / __// __ \ / ___// _ \
 / ____// /   /  __// /__ / /_/ // // /__ / /_/ // // /_/ // /_ / // /_/ // / / /   ___/ // /_ / /_/ // /   /  __/
/_/    /_/    \___/ \___/ \__,_//_/ \___/ \__,_//_/ \__,_/ \__//_/ \____//_/ /_/   /____/ \__/ \____//_/    \___/                                                     
"""


def echo_intro(click_iface, mode):
    width = 120
    mode = (
        "Hybrid [Local + Cloud]"
        if mode == "cache-only"
        else mode.replace("-", " ").capitalize()
    )
    click_iface.echo(f"{title:^{width}}", fg="red", bold=True)
    mode_text = f"⚙️ Fetch Mode: {mode}"
    version_text = "🏷️ [Version 0.1.0] 🏷️"
    click_iface.echo(f"{version_text:^{width}}", fg="cyan", bold=True)
    click_iface.echo(f"{mode_text:^{width}}", fg="cyan", bold=True)
    click_iface.echo("")


def echo_exceptions(message, click_iface, exit=False, bg="red", fg="white"):
    click_iface.echo(f"{message}", fg=fg, bg=bg, bold=False)
    if exit:
        echo_sys_exited(click_iface)
        sys.exit(1)


def echo_redis_job_submitted(click_iface, message):
    click_iface.echo(
        f"{log_prefix()}Job submitted to Redis caching system. {message}",
        fg="blue",
        bold=False,
    )


def echo_redis_fetched_missed(click_iface, size_res, size_missing):
    click_iface.echo(
        f"{log_prefix()}Results size of {size_res} fetched from local cache. Missing inputs: {size_missing}",
        fg="blue",
        bold=False,
    )


def echo_redis_local_completed(click_iface):
    click_iface.echo(
        f"{log_prefix()}All inputs has results in local caches and won't processed further query.",
        fg="blue",
        bold=False,
    )
    echo_sys_exited(click_iface)


def echo_uploading_inputs(click_iface):
    click_iface.echo(
        f"{log_prefix()}Uploading input data to S3 bucket", fg="blue", bold=False
    )


def echo_local_fetched_cache_szie(click_iface, cache_size, none_count):
    if cache_size != 0:
        click_iface.echo(
            f"{log_prefix()}Found cache size of {cache_size} from local Redis cache. Post processing started!",
            fg="blue",
            bold=False,
        )
    else:
        echo_local_only_empty_cache(click_iface)
        echo_redis_null_output(click_iface)
    if none_count != 0:
        echo_redis_null_output(click_iface)


def echo_upload_complete(click_iface):
    click_iface.echo(f"{log_prefix()}Upload completed.", fg="blue", bold=False)


def echo_submitting_job(click_iface, model_id):
    click_iface.echo(
        f"{log_prefix()}Submitting a task for model[{model_id}] precalculation",
        fg="blue",
        bold=False,
    )


def echo_sys_exited(click_iface):
    click_iface.echo(
        f"{log_prefix()}System exits with exit code 1",
        fg="blue",
        bold=False,
    )


def echo_redis_null_output(click_iface):
    click_iface.echo(
        f"{log_prefix()}Beware that output file may contain None or empty values!",
        fg="yellow",
        bold=True,
    )


def echo_job_submitted(click_iface, job_id: str):
    click_iface.echo(
        f"{log_prefix()}Job submitted with a job id: {job_id}", fg="blue", bold=False
    )


def echo_local_only_empty_cache(click_iface):
    click_iface.echo(
        f"{log_prefix()}Ooopps! No cached results found in local Redis system for the model! Try cloud caching!",
        fg="white",
        blink=False,
        bold=True,
        bg="yellow",
    )


def echo_local_sample_warning_(
    click_iface, n: int, cache_size: int, output_path: str = None
) -> bool:
    prompt = "Do you want to continue to cloud for fetching?"

    if cache_size == 0:
        echo_local_only_empty_cache(click_iface)
        return click.confirm(prompt)

    d = n - cache_size
    if cache_size >= n:
        bg = "cyan"
        message = "This is more or equal to the sample size you requested!"
    else:
        bg = "yellow"
        message = f"This is less than the sample size you requested by {d}!"

    click_iface.echo(
        f"{log_prefix()}Cache size of {cache_size} fetched from local Redis caching. {message}",
        fg="white",
        blink=False,
        bold=True,
        bg=bg,
    )
    if cache_size >= n:
        echo_redis_file_saved(click_iface, output_path)
        return False
    return click.confirm(prompt)


def echo_small_sample_warning(click_iface, n: int):
    if n <= 50000:
        promompt = "Do you want to continue?"
        click_iface.echo(
            f"{log_prefix()}Sample size of less than 50,000 [{n}] is not recommended for fetching precalculation!",
            fg="white",
            blink=True,
            bg="red",
        )
        return click.confirm(promompt)


def echo_status(click_iface, status: str):
    if status == JobStatus.PENDING:
        click_iface.echo(
            f"{log_prefix()}Submitted job status: {status}", fg="yellow", bold=False
        )
    elif status == JobStatus.FAILED:
        click_iface.echo(
            f"{log_prefix()}Submitted job status: {status}", fg="red", bold=False
        )
    elif status == JobStatus.RUNNING:
        click_iface.echo(
            f"{log_prefix()}Submitted job status: {status}", fg="cyan", bold=False
        )
    else:
        click_iface.echo(
            f"{log_prefix()}Submitted job status: {status}", fg="green", bold=False
        )


def echo_job_succeeded(click_iface):
    click_iface.echo(
        f"{log_prefix()}Precalculation successfully fetched", fg="blue", bold=False
    )


def echo_found_shards(click_iface, count: int):
    click_iface.echo(
        f"{log_prefix()}Fetched {count} .gz compressed precalcultaion files.",
        fg="blue",
        bold=False,
    )


def echo_merged_saved(click_iface, output_path: Path):
    click_iface.echo(
        f"{log_prefix()}Merged CSV saved to: {output_path}", fg="blue", bold=False
    )


def echo_redis_file_saved(click_iface, output_path: Path):
    click_iface.echo(
        f"{log_prefix()}Output file is saved to: {output_path}", fg="blue", bold=False
    )


def merge_csvs_stdlib(files, output):
    with open(output, "w", newline="") as fout:
        writer = None
        for path in files:
            with open(path, "r", newline="") as fin:
                reader = csv.reader(fin)
                header = next(reader)
                if writer is None:
                    writer = csv.writer(fout)
                    writer.writerow(header)
                for row in reader:
                    writer.writerow(row)


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
            writer.writerow([rec["key"], rec["input"]])
        tmp.flush()
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
        shards,
        header,
        smiles_list,
        output_path,
        downloader,
        click_iface,
        local_cache_csv_path,
    ):
        """
        Download, decompress, and concatenate multiple CSV shards into a single file.
        If `smiles_list` is non-empty, reorder the rows to match its order;
        otherwise, just merge them sequentially.

        """
        gz_shards = [
            s
            for s in shards
            if s.get("key", "").endswith(".gz") or s.get("url", "").endswith(".gz")
        ]
        total_size = 0
        for s in gz_shards:
            size = s.get("size", 0)
            try:
                hdrs = downloader.head(s["url"])
                size = int(hdrs.get("Content-Length", size))
            except Exception:
                pass
            total_size += size

        click_iface.echo(
            f"⬇ Total download size: {total_size / 1024:.1f} KB", fg="blue"
        )

        if not smiles_list:
            click_iface.echo(
                "⬇ No reorder list provided; merging sequentially", fg="blue"
            )
            with open(output_path, "w", newline="", encoding="utf-8") as out_f:
                writer = csv.writer(out_f)
                writer.writerow(header)
                pbar = click_iface.progress_bar(
                    total=total_size,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc="⬇ downloading shards",
                    colour="blue",
                    ncols=80,
                )
                first_shard = True
                for shard in gz_shards:
                    buf = io.BytesIO()
                    for chunk in downloader.download_stream(shard["url"]):
                        buf.write(chunk)
                        pbar.update(len(chunk))
                    text = gzip.decompress(buf.getvalue()).decode(
                        "utf-8", errors="replace"
                    )
                    reader = csv.reader(io.StringIO(text))
                    if first_shard:
                        first_shard = False
                    for row in reader:
                        writer.writerow(row)
                pbar.close()

            if (
                os.path.exists(local_cache_csv_path)
                and os.path.getsize(local_cache_csv_path) > 0
            ):
                # fmt: off
                with open(local_cache_csv_path, "r", encoding="utf-8") as cache_f, open(output_path, "a", newline="", encoding="utf-8") as out_f:
                    next(cache_f, None)
                    for line in cache_f:
                        out_f.write(line)
                try:
                    os.remove(local_cache_csv_path)
                except OSError:
                    pass
                # fmt: off

            click_iface.echo(
                f"✅ Sequential merge written to {output_path}", fg="green"
            )
            return

        if "input" not in header:
            raise KeyError("Header must contain 'input' to reorder")

        col_idx = header.index("input")

        normalized = []
        for item in smiles_list:
            if isinstance(item, str):
                normalized.append(item)
            elif isinstance(item, dict) and "input" in item:
                normalized.append(item["input"])
            else:
                raise TypeError(f"Cannot extract 'input' from {item}")

        click_iface.echo("⬇ Building lookup and reordering by 'input'", fg="blue")
        lookup = {}
        pbar = click_iface.progress_bar(
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc="⬇ downloading shards",
            colour="blue",
            ncols=80,
        )
        first_shard = True
        for shard in gz_shards:
            buf = io.BytesIO()
            for chunk in downloader.download_stream(shard["url"]):
                buf.write(chunk)
                pbar.update(len(chunk))
            text = gzip.decompress(buf.getvalue()).decode("utf-8", errors="replace")
            reader = csv.reader(io.StringIO(text))
            if first_shard:
                first_shard = False
            for row in reader:
                row_len = len(row) // 2
                row = row[:row_len] if row_len > 1 else row
                lookup[row[col_idx]] = row
        pbar.close()

        with open(output_path, "w", newline="", encoding="utf-8") as out_f:
            writer = csv.writer(out_f)
            writer.writerow(header)
            for inp in normalized:
                row = lookup.get(inp)
                if row:
                    writer.writerow(row)
                else:
                    empty_row = [""] * len(header)
                    empty_row[col_idx] = inp
                    writer.writerow(empty_row)

        if (
            os.path.exists(local_cache_csv_path)
            and os.path.getsize(local_cache_csv_path) > 0
        ):
            with open(local_cache_csv_path, "r", encoding="utf-8") as cache_f:
                cache_reader = csv.reader(cache_f)
                next(cache_reader, None)
                cache_rows = list(cache_reader)

            with open(output_path, "r", newline="", encoding="utf-8") as merged_f:
                merged_reader = csv.reader(merged_f)
                merged_header = next(merged_reader)
                merged_rows = list(merged_reader)

            if "input" in merged_header:
                idx = merged_header.index("input")
            else:
                idx = None

            def is_full(row):
                return all(cell.strip() for cell in row)

            final_cache = []
            for cro in cache_rows:
                match = None
                for mro in merged_rows:
                    if idx is not None and cro[idx] == mro[idx]:
                        match = mro
                        break
                    elif idx is None and cro == mro:
                        match = mro
                        break

                if match and not is_full(cro) and is_full(match):
                    final_cache.append(match)
                else:
                    final_cache.append(cro)

            with open(output_path, "w", newline="", encoding="utf-8") as wf:
                writer = csv.writer(wf)
                writer.writerow(merged_header)

                seen = set()
                for row in final_cache:
                    writer.writerow(row)
                    key = row[idx] if idx is not None else tuple(row)
                    seen.add(key)

                for mro in merged_rows:
                    key = mro[idx] if idx is not None else tuple(mro)
                    if key not in seen:
                        writer.writerow(mro)
            try:
                os.remove(local_cache_csv_path)
            except OSError:
                pass

        click_iface.echo(f"✅ Final ordered CSV written to {output_path}", fg="green")


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
    CACHE_ONLY = "cache-only"
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

    @classmethod
    def is_precalculation_enabled(cls, option):
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
        return (
            option == cls.CLOUD_ONLY
            or option == cls.LOCAL_ONLY
            or option == cls.CACHE_ONLY
        )


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
