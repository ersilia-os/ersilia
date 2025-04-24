import csv
import gzip
import io
import tempfile
import time
import uuid
from pathlib import Path

import click
import requests
from tqdm import tqdm

from ersilia.core.base import ErsiliaBase
from ersilia.default import API_BASE, INFERENCE_STORE_API_URL
from ersilia.io.input import GenericInputAdapter
from ersilia.io.output import GenericOutputAdapter
from ersilia.store.utils import OutputSource, delete_file_upon_upload


class InferenceStoreApi(ErsiliaBase):
    """
    Synchronous client for Inference Store API. Uploads inputs, submits/polls a job,
    fetches shard URLs, and downloads & merges all CSV shards into one UTF-8 CSV with headers,
    while displaying click-styled info and a tqdm progress bar.
    """

    def __init__(
        self,
        model_id: str,
        output: str = None,
        n_samples: int = -1,
        output_source: str = "cloud-cache-only",
    ):
        super().__init__()
        self.model_id = model_id
        self.n_samples = n_samples
        self.output_source = output_source
        self.request_id = None
        # Pre-fetch output schema
        cols = GenericOutputAdapter(model_id=model_id)._fetch_schema_from_github()[0]
        self.header = ["key", "input"] + cols
        self.output_path = Path(output) if output else Path(f"{model_id}_precalc.csv")
        self.input_adapter = GenericInputAdapter(model_id=model_id)

    def _gen_request_id(self):
        self.request_id = str(uuid.uuid4())

    def _upload_inputs(self, inputs):
        click.echo(click.style("Uploading input CSV to S3...", fg="green", bold=True))
        # 1) get presigned URL
        resp = requests.get(
            f"{INFERENCE_STORE_API_URL}/upload-destination",
            params={"modelid": self.model_id, "requestid": self.request_id},
            timeout=60,
        )
        resp.raise_for_status()
        pres = resp.json()

        # 2) write temp CSV
        tmp = tempfile.NamedTemporaryFile(
            "w+", delete=False, newline="", encoding="utf-8"
        )
        writer = csv.writer(tmp)
        for rec in self.input_adapter.adapt_one_by_one(inputs):
            writer.writerow([rec["input"]])
        tmp.flush()
        tmp.close()

        # 3) upload file
        with open(tmp.name, "rb") as fbin:
            files = {"file": (tmp.name, fbin)}
            up = requests.post(
                pres["url"], data=pres.get("fields", {}), files=files, timeout=300
            )
        delete_file_upon_upload(up.status_code, tmp.name)
        up.raise_for_status()
        click.echo(click.style("Upload complete.", fg="green", bold=True))

    def _submit_job(self) -> str:
        click.echo(click.style("Submitting Athena UNLOAD job...", fg="cyan", bold=True))
        payload = {
            "requestid": self.request_id,
            "modelid": self.model_id,
            "fetchtype": "all",
            "nsamples": self.n_samples,
        }
        r = requests.post(f"{API_BASE}/submit", json=payload, timeout=10)
        r.raise_for_status()
        job_id = r.json()["jobId"]
        click.echo(click.style(f"Job submitted: {job_id}", fg="cyan"))
        return job_id

    def _poll_status(self, job_id: str, interval: float = 5.0, timeout: float = 900.0):
        click.echo(click.style("Polling job status...", fg="yellow", bold=True))
        start = time.time()
        while time.time() - start < timeout:
            r = requests.get(f"{API_BASE}/status", params={"jobId": job_id}, timeout=10)
            r.raise_for_status()
            st = r.json().get("status")
            click.echo(click.style(f"Status: {st}", fg="yellow"))
            if st == "SUCCEEDED":
                click.echo(click.style("Job succeeded!", fg="green", bold=True))
                return
            if st in ("FAILED", "ERROR"):
                msg = r.json().get("errorMessage", "Unknown error")
                raise RuntimeError(f"Job failed: {msg}")
            time.sleep(interval)
        raise RuntimeError("Job polling timed out")

    def _fetch_shards(self, job_id: str) -> list[dict]:
        click.echo(click.style("Fetching result shard URLs...", fg="blue", bold=True))
        r = requests.get(f"{API_BASE}/result", params={"jobId": job_id}, timeout=10)
        r.raise_for_status()
        files = r.json().get("files", [])
        if not files:
            raise RuntimeError("No shards returned")
        click.echo(click.style(f"Found {len(files)} shards.", fg="blue"))
        return files

    def get_precalculations(self, inputs: list) -> str:
        """
        Full flow: upload inputs (if needed), submit job, poll, calculate total size,
        show progress bar while downloading & merging shards into one CSV.
        Returns the path to combined UTF-8 CSV with headers.
        """
        self._gen_request_id()
        if self.output_source == OutputSource.CLOUD:
            self._upload_inputs(inputs)

        job_id = self._submit_job()
        self._poll_status(job_id)
        shards = self._fetch_shards(job_id)

        # Calculate total download size via HEAD
        total_size = 0
        click.echo(click.style("Calculating total download size...", fg="magenta"))
        for shard in shards:
            h = requests.head(shard["url"], timeout=10)
            size = int(h.headers.get("Content-Length", 0))
            total_size += size
        click.echo(
            click.style(f"Total download: {total_size/1024:.1f} KB", fg="magenta")
        )

        # Download & merge with progress bar
        click.echo(
            click.style("Downloading and merging shards...", fg="green", bold=True)
        )
        with open(self.output_path, "w", encoding="utf-8", newline="") as out_f:
            out_f.write(",".join(self.header) + "\n")
            first = True
            pbar = tqdm(
                total=total_size,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc="⬇ CSV",
                colour="blue",
                ncols=70,
                bar_format="{desc}: [{bar}] {percentage:3.0f}% {n_fmt}/{total_fmt} • {rate_fmt} • ETA: {remaining}",
            )
            for shard in shards:
                url = shard.get("url")
                key = shard.get("key", "")
                r = requests.get(url, stream=True, timeout=300)
                r.raise_for_status()
                buf = io.BytesIO()
                for chunk in r.iter_content(64 * 1024):
                    buf.write(chunk)
                    pbar.update(len(chunk))
                data = buf.getvalue()
                if key.endswith(".gz") or url.endswith(".gz"):
                    data = gzip.decompress(data)
                text = data.decode("utf-8", errors="replace")
                for i, line in enumerate(text.splitlines(keepends=True)):
                    if not first and i == 0 and line.startswith(self.header[0]):
                        continue
                    out_f.write(line)
                first = False
            pbar.close()
        click.echo(
            click.style(
                f"Merged CSV saved to: {self.output_path}", fg="green", bold=True
            )
        )
        return str(self.output_path)
