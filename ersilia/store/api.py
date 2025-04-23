import asyncio
import csv
import tempfile
import time
import uuid
from pathlib import Path

import aiohttp
import click
import requests
from aiohttp import ClientConnectorError, ClientResponseError, ClientTimeout
from tqdm.asyncio import tqdm

from ersilia.core.base import ErsiliaBase
from ersilia.default import API_BASE, CLOUD_CACHE_CHUNK, INFERENCE_STORE_API_URL
from ersilia.io.input import GenericInputAdapter
from ersilia.store.utils import (
    OutputSource,
    delete_file_upon_upload,
)


class InferenceStoreApi(ErsiliaBase):
    """
    A class to interact with the Inference Store API.

    This class provides methods to upload input data to the inference store,
    request predictions, and retrieve the results. It handles the creation
    of unique request IDs, obtaining presigned URLs for file uploads, and
    managing the input data format.

    Parameters
    ----------
    model_id : str
        The ID of the model for which the inference is being performed.

    Examples
    --------
    .. code-block:: python

        api = InferenceStoreApi(model_id="eosxxxx")
        inputs = ["CCO", "CCN"]
        api.get_precalculations(inputs)

    Mechanism
    ---------
    1. Initialize the class with a model ID.
    2. Generate a unique request ID for each inference request.
    3. Obtain a presigned URL from the inference store for uploading input data.
    4. Adapt the input data to the required format and upload it using the presigned URL.
    5. Request the predictions from the inference store using the request ID.
    6. Retrieve and return the URL of the output predictions.
    """

    def __init__(
        self,
        model_id: str,
        output: str,
        output_source: str = "cloud-cache-only",
        n_samples: int = -1,
    ):
        ErsiliaBase.__init__(self)
        self.n_samples = n_samples
        self.model_id = model_id
        self.output_source = output_source
        self.output_path = Path(f"{self.model_id}_cloud_cache.csv")
        self.output = output if output is not None else self.output_path
        self.request_id = None
        self.input_adapter = GenericInputAdapter(model_id=self.model_id)
        self._connector = None
        self._session = None

    async def _get_session(
        self, connector_kwargs: dict | None = None, session_kwargs: dict | None = None
    ) -> aiohttp.ClientSession:
        default_connector = {"limit": 50, "limit_per_host": 10}
        default_session = {}

        connector_kwargs = {**default_connector, **(connector_kwargs or {})}
        session_kwargs = {**default_session, **(session_kwargs or {})}

        if self._session is None or self._session.closed:
            self._connector = aiohttp.TCPConnector(**connector_kwargs)
            self._session = aiohttp.ClientSession(
                connector=self._connector, **session_kwargs
            )
        return self._session

    def _generate_request_id(self) -> str:
        self.request_id = str(uuid.uuid4())
        self.logger.debug(f"Generating the request id: {self.request_id}")

    def _get_presigned_url(self) -> str:
        response = requests.get(
            INFERENCE_STORE_API_URL + "/upload-destination",
            params={"modelid": self.model_id, "requestid": self.request_id},
        )

        presigned_url_response = response.json()
        self.logger.debug(
            f"Generating the pre-assigned url: {presigned_url_response['url'][:12]}..."
        )
        return presigned_url_response

    def _post_inputs(self, inputs, presigned_url_response) -> str:
        adapted_input_generator = self.input_adapter.adapt_one_by_one(inputs)
        self.logger.info(
            "Converting inputs to smiles list for posting to inference store!"
        )
        smiles_list = [val["input"] for val in adapted_input_generator]
        self.logger.debug(
            f"Smiles with the size of {len(smiles_list)} has been processed!"
        )
        # TODO: avoid creating a file
        ersilia_input_file = "ersilia_cloud_inputs.csv"
        with open(ersilia_input_file, "w", newline="") as input_file:
            writer = csv.writer(input_file)
            for smiles in smiles_list:
                writer.writerow([smiles])

        with open(ersilia_input_file, "rb") as input_file:
            files = {"file": (ersilia_input_file, input_file)}
            presigned_url = presigned_url_response.get("url")
            presigned_url_data = presigned_url_response.get("fields")
            self.logger.info("Posting input to inference and fetching results!")

            http_response = requests.post(
                presigned_url, data=presigned_url_data, files=files, timeout=10000
            )

        self.logger.debug("Deleting file upon uploading!")
        delete_file_upon_upload(http_response.status_code, ersilia_input_file)
        return http_response

    def _submit_job(self) -> str:
        url = API_BASE + "/submit"
        payload = {
            "requestid": self.request_id,
            "modelid": self.model_id,
            "fetchtype": "all",
            "nsamples": self.n_samples,
        }

        response = requests.post(url, json=payload, timeout=10)
        response.raise_for_status()
        data = response.json()
        return data["jobId"]

    def _get_presigned_url_from_s3(
        self, job_id: str, poll_interval: float = 20.0, timeout: float = 1440000
    ) -> str:
        status_url = API_BASE + "/status"
        result_url = API_BASE + "/result"
        start_time = time.time()

        while True:
            elapsed = time.time() - start_time
            if elapsed > timeout:
                raise Exception("Timed out waiting for job completion")

            params = {"jobId": job_id}
            status_response = requests.get(status_url, params=params, timeout=10)
            status_response.raise_for_status()
            status_data = status_response.json()

            current_status = status_data.get("status")
            click.echo(
                click.style(
                    f"Job {job_id} current status: {current_status}",
                    fg="green",
                    bold=True,
                )
            )

            if current_status == "SUCCEEDED":
                break
            elif current_status in ["FAILED", "ERROR"]:
                error_message = status_data.get("errorMessage", "Unknown error")
                raise Exception(
                    f"Job {job_id} failed with status '{current_status}': {error_message}"
                )

            time.sleep(poll_interval)

        result_response = requests.get(result_url, params={"jobId": job_id}, timeout=10)
        result_response.raise_for_status()
        result_data = result_response.json()
        presigned_url = result_data.get("url")
        size = result_data.get("size")
        if not presigned_url:
            raise Exception("No URL returned in result")
        return presigned_url, size

    def _get_outputs(self) -> str:
        self.logger.info("Getting result output from a preassigned url")
        job_id = self._submit_job()
        self.logger.info(
            f"Job id [jobId] fetched from a job submission endpoint: {job_id}"
        )
        self.logger.debug(
            "Waiting the job process to get preassigned url for prediction file!"
        )
        url, size = self._get_presigned_url_from_s3(job_id)
        self.logger.info(f"Successfully fetched the prediction file url:{url[:17]}")
        return url, size

    def _get_precalculation(self, inputs: list) -> str:
        self._generate_request_id()
        if self.output_source == OutputSource.CLOUD:
            presigned_url = self._get_presigned_url()
            response = self._post_inputs(inputs, presigned_url)
            if response.status_code == 204:
                self.logger.info("File uploaded successfully")
            else:
                return f"Failed to upload file: {response.status_code} error ({response.text})"
        try:
            st = time.perf_counter()
            output_presigned_url, size = self._get_outputs()
        except RuntimeError as e:
            return f"Error retrieving outputs: {str(e)}"
        if not output_presigned_url:
            self.logger.error("No outputs found in store!")
            return None
        et = time.perf_counter()
        self.logger.info(
            f"Response fetched from the inference store in: {et - st:.4f} seconds!"
        )
        return output_presigned_url, size

    async def _fetch_csv_async(self, chunk_inputs: list) -> str:
        raw, size = self._get_precalculation(chunk_inputs)
        self.logger.info(f"Size: {size}")
        if isinstance(raw, dict):
            url = raw.get("url") or raw.get("URL") or ""
            size = size.get("size") or size.get("SIZE") or ""
        else:
            url = raw if isinstance(raw, str) else ""

        if not url or url.startswith(("Failed", "Error")):
            self.logger.error(f"Invalid precalc URL: {raw!r}")
            return ""

        timeout = ClientTimeout(total=None, sock_connect=30, sock_read=30)
        session = await self._get_session(
            connector_kwargs={"limit_per_host": 4}, session_kwargs={"timeout": timeout}
        )
        self.logger.info("Starting download of one chunk…")

        tmp_path = Path(tempfile.gettempdir()) / f"{uuid.uuid4()}.csv"
        pbar = tqdm(
            total=size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc="⬇ CSV",
            colour="blue",
            ncols=70,
            bar_format="{desc}: [{bar}] {percentage:3.0f}% {n_fmt}/{total_fmt} • {rate_fmt} • ETA: {remaining}",
        )

        try:
            async with session.get(url) as resp:
                resp.raise_for_status()
                with open(tmp_path, "wb") as f:
                    async for block in resp.content.iter_chunked(64 * 1024):
                        f.write(block)
                        pbar.update(len(block))

        except (ClientConnectorError, ClientResponseError, asyncio.TimeoutError) as e:
            self.logger.error(f"Chunk download failed: {e!r}")
            if tmp_path.exists():
                tmp_path.unlink()
            return ""
        finally:
            pbar.close()
            await session.close()
            if self._connector:
                self._connector.close()
            self._session = self._connector = None

        return str(tmp_path)

    def _fetch_csv_sync(self, chunk_inputs: list) -> str:
        try:
            return asyncio.run(self._fetch_csv_async(chunk_inputs))
        except Exception as e:
            self.logger.error(f"Error in synchronous wrapper: {e}")
            return ""

    def _echo(self, text, fg="green"):
        click.echo(
            click.style(
                text,
                fg=fg,
                bold=True,
            )
        )

    def _get_precalculations(self, inputs: list) -> str:
        self.logger.warning("Starting full CSV precalc…")
        start = time.perf_counter()

        temp_files = []
        if inputs:
            for i in range(0, len(inputs), CLOUD_CACHE_CHUNK):
                chunk = inputs[i : i + CLOUD_CACHE_CHUNK]
                self._echo(f"Downloading chunk #{i // CLOUD_CACHE_CHUNK + 1}")
                path = self._fetch_csv_sync(chunk)
                if path:
                    temp_files.append(path)
        else:
            self._echo("Downloading all precalculations", fg="cyan")
            path = self._fetch_csv_sync(None)
            if path:
                temp_files.append(path)

        self._echo(f"All chunks downloaded in {time.perf_counter() - start:.2f}s")
        if not temp_files:
            self.logger.error("No chunks → aborting")
            return "No data to write."

        writer = None

        with open(self.output, "w", newline="") as out_f:
            for fp in temp_files:
                with open(fp, newline="") as in_f:
                    reader = csv.DictReader(in_f)

                    if writer is None:
                        writer = csv.DictWriter(out_f, fieldnames=reader.fieldnames)
                        writer.writeheader()

                    for row in reader:
                        writer.writerow(row)

        self.logger.info(f"Combined CSV written to {self.output}")
        return self.output

    def get_precalculations(self, inputs: list = None) -> str:
        """
        To fetch precalculation
        """
        click.echo(
            click.style(
                f"Fetching precalculation for model: {self.model_id} started. Please wait patiently!",
                bg="yellow",
                fg="white",
            )
        )
        res = self._get_precalculations(inputs)
        return res
