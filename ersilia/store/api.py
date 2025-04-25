import time
import uuid
from pathlib import Path

from ersilia.core.base import ErsiliaBase
from ersilia.default import API_BASE, INFERENCE_STORE_API_URL
from ersilia.io.input import GenericInputAdapter
from ersilia.io.output import GenericOutputAdapter
from ersilia.store.utils import (
    ApiClient,
    ClickInterface,
    FileManager,
    JobStatus,
    OutputSource,
    echo_found_shards,
    echo_job_submitted,
    echo_job_succeeded,
    echo_merged_saved,
    echo_status,
    echo_submitting_job,
    echo_upload_complete,
    echo_uploading_inputs,
)


class InferenceStoreApi(ErsiliaBase):
    """
    Synchronous client for the Inference Store API orchestrating end-to-end batch cloud precalculation fetcher.

    Parameters
    ----------
    model_id : str
        Identifier of the model to query.
    output : str, optional
        Path to save the merged CSV (default is "<model_id>_precalc.csv").
    n_samples : int, optional
        Number of samples to retrieve (default is -1 for all available).
    output_source : OutputSource, optional
        Source of input data (default is OutputSource.CLOUD_ONLY).
    click_iface : ClickInterface, optional
        CLI helper instance for styled output (default creates new one).
    api_client : ApiClient, optional
        HTTP client for API interactions (default creates new one).
    file_manager : FileManager, optional
        File operations helper (default creates new one).
    """

    def __init__(
        self,
        model_id: str,
        output: str = None,
        n_samples: int = -1,
        output_source: str = OutputSource.CLOUD_ONLY,
        click_iface: ClickInterface = None,
        api_client: ApiClient = None,
        file_manager: FileManager = None,
    ):
        super().__init__()
        self.model_id = model_id
        self.n_samples = n_samples
        self.output_source = output_source
        self.request_id = None
        cols = GenericOutputAdapter(model_id=model_id)._fetch_schema_from_github()[0]
        self.header = ["key", "input"] + cols
        self.output_path = Path(output) if output else Path(f"{model_id}_precalc.csv")
        self.input_adapter = GenericInputAdapter(model_id=model_id)
        self.click = click_iface or ClickInterface()
        self.api = api_client or ApiClient()
        self.files = file_manager or FileManager()

    def get_precalculations(self, inputs: list) -> str:
        """
        Execute full inference workflow: upload inputs, submit job, poll status,
        retrieve result shards, and merge them into a single CSV.

        Parameters
        ----------
        inputs : list
            List of input records to process.

        Returns
        -------
        str
            File path to the merged CSV containing inference results.

        Raises
        ------
        RuntimeError
            If the job fails, no shards are returned, or polling times out.
        """
        self.request_id = str(uuid.uuid4())
        if self.output_source == OutputSource.CLOUD:
            echo_uploading_inputs(self.click)
            pres = self.api.get_json(
                f"{INFERENCE_STORE_API_URL}/upload-destination",
                params={"modelid": self.model_id, "requestid": self.request_id},
            )
            tmp_path = self.files.create_temp_csv(inputs, self.input_adapter)
            self.files.upload_to_s3(pres, tmp_path)
            echo_upload_complete(self.click)

        echo_submitting_job(self.click)
        payload = {
            "requestid": self.request_id,
            "modelid": self.model_id,
            "fetchtype": "all",
            "nsamples": self.n_samples,
        }
        job_id = self.api.post_json(f"{API_BASE}/submit", json=payload)["jobId"]
        echo_job_submitted(self.click, job_id)

        start = time.time()
        while True:
            status = self.api.get_json(f"{API_BASE}/status", params={"jobId": job_id})[
                "status"
            ]
            echo_status(self.click, status)
            if status == JobStatus.SUCCEEDED:
                echo_job_succeeded(self.click)
                break
            if status == JobStatus.FAILED:
                error = self.api.get_json(
                    f"{API_BASE}/status", params={"jobId": job_id}
                ).get("errorMessage", "Unknown error")
                raise RuntimeError(f"Job failed: {error}")
            if time.time() - start > 3600:
                raise RuntimeError("Job polling timed out")
            time.sleep(5)

        shards = self.api.get_json(f"{API_BASE}/result", params={"jobId": job_id})[
            "files"
        ]
        if not shards:
            raise RuntimeError("No shards returned")
        echo_found_shards(self.click, len(shards))

        self.files.merge_shards(
            shards, self.header, self.output_path, self.api, self.click
        )

        echo_merged_saved(self.click, self.output_path)
        return str(self.output_path)
