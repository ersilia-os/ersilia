import csv
import gzip
import tempfile
import time
import uuid
from pathlib import Path

import requests

from ersilia.core.base import ErsiliaBase
from ersilia.default import API_BASE, INFERENCE_STORE_API_URL
from ersilia.io.input import GenericInputAdapter
from ersilia.io.output import GenericOutputAdapter
from ersilia.store.utils import OutputSource, delete_file_upon_upload


class InferenceStoreApi(ErsiliaBase):
    """
    Synchronous client for Inference Store API. Uploads inputs, submits/polls a job,
    fetches shard URLs, downloads shards into memory, decompresses, and merges
    into a single UTF-8 CSV with proper headers.
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

    def _submit_job(self) -> str:
        payload = {
            "requestid": self.request_id,
            "modelid": self.model_id,
            "fetchtype": "all",
            "nsamples": self.n_samples,
        }
        r = requests.post(f"{API_BASE}/submit", json=payload, timeout=10)
        r.raise_for_status()
        return r.json()["jobId"]

    def _poll_status(self, job_id: str, interval: float = 5.0, timeout: float = 900.0):
        start = time.time()
        while time.time() - start < timeout:
            r = requests.get(f"{API_BASE}/status", params={"jobId": job_id}, timeout=10)
            r.raise_for_status()
            st = r.json().get("status")
            if st == "SUCCEEDED":
                return
            if st in ("FAILED", "ERROR"):
                raise RuntimeError(f"Job failed: {r.json().get('errorMessage')}")
            time.sleep(interval)
        raise RuntimeError("Job timed out")

    def _fetch_shards(self, job_id: str) -> list[dict]:
        r = requests.get(f"{API_BASE}/result", params={"jobId": job_id}, timeout=10)
        r.raise_for_status()
        files = r.json().get("files", [])
        if not files:
            raise RuntimeError("No shards returned")
        return files

    def get_precalculations(self, inputs: list) -> str:
        """
        Full flow: upload inputs (if needed), submit job, poll, download all shards in-memory,
        decompress, merge with headers, write UTF-8 CSV file, return its path.
        """
        self._gen_request_id()
        if self.output_source == OutputSource.CLOUD:
            self._upload_inputs(inputs)

        job_id = self._submit_job()
        self._poll_status(job_id)
        shards = self._fetch_shards(job_id)

        # Combine shards
        lines_written = 0
        with open(self.output_path, "w", encoding="utf-8", newline="") as out_f:
            # write header
            out_f.write(",".join(self.header) + "\n")
            lines_written += 1

            for shard in shards:
                url = shard.get("url")
                key = shard.get("key", "")
                # download entire shard
                r = requests.get(url, timeout=300)
                r.raise_for_status()
                data = r.content  # bytes

                # decompress if gzipped
                if key.endswith(".gz") or shard.get("url", "").endswith(".gz"):
                    data = gzip.decompress(data)

                # decode to UTF-8, replace errors
                text = data.decode("utf-8", errors="replace")
                # split into lines keeping newline
                for i, line in enumerate(text.splitlines(keepends=True)):
                    # skip header in subsequent shards
                    if (
                        lines_written > 1
                        and i == 0
                        and line.rstrip().split(",") == self.header
                    ):
                        continue
                    out_f.write(line)
                    lines_written += 1

        return str(self.output_path)
