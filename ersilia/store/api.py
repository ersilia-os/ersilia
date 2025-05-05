import json
import sys
import time
import uuid
from pathlib import Path

from ersilia.core.base import ErsiliaBase
from ersilia.default import API_BASE, DEFAULT_API_NAME, EOS_TMP, INFERENCE_STORE_API_URL
from ersilia.io.input import GenericInputAdapter
from ersilia.io.output import GenericOutputAdapter
from ersilia.store.dump import DumpLocalCache
from ersilia.store.utils import (
    ApiClient,
    ClickInterface,
    FileManager,
    JobStatus,
    OutputSource,
    echo_found_shards,
    echo_intro,
    echo_job_submitted,
    echo_job_succeeded,
    echo_local_fetched_cache_szie,
    echo_local_only_empty_cache,
    echo_local_sample_warning,
    echo_merged_saved,
    echo_redis_fetched_missed,
    echo_redis_file_saved,
    echo_redis_job_submitted,
    echo_redis_local_completed,
    echo_small_sample_warning,
    echo_status,
    echo_submitting_job,
    echo_sys_exited,
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
        self.dump_local = DumpLocalCache()
        self.fetch_type = "all"
        self.generic_output_adapter = GenericOutputAdapter(model_id=model_id)
        self.schema = self.generic_output_adapter._fetch_schema_from_github()
        assert self.schema is not None, "Model schema can not be fetched from github."
        self.cols = self.schema[0]
        self.dtype = self.schema[1]
        self.header = ["key", "input"] + self.cols
        self.output_path = Path(output) if output else Path(f"{model_id}_precalc.csv")
        self.input_adapter = GenericInputAdapter(model_id=model_id)
        self.click = click_iface or ClickInterface()
        self.api = api_client or ApiClient()
        self.files = file_manager or FileManager()
        self.local_cache_csv_path = f"{EOS_TMP}/local.csv"

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
        echo_intro(self.click, self.output_source)
        if self.output_source == OutputSource.LOCAL_ONLY:
            return self._handle_strict_local(inputs)

        shards = self._submit_and_get_shards(inputs)
        return self._merge_shards(shards)

    def _handle_strict_local(self, inputs: list) -> str:
        self.dump_local.init_redis()
        if inputs:
            echo_redis_job_submitted(self.click, f"Input size: {len(inputs)}")
            results = self.dump_local.get_cached(
                self.model_id, inputs, self.dtype, cols=self.cols
            )
            cache_size = abs(len(results[0]) - len(results[1]))
            echo_local_fetched_cache_szie(self.click, cache_size)
        else:
            ns = self.n_samples if self.n_samples != -1 else "all"
            echo_redis_job_submitted(self.click, f"Sample size: {ns}")
            results, inputs = self.dump_local.fetch_all_cached(
                self.model_id, self.dtype, cols=self.cols
            )
            if not results:
                echo_local_only_empty_cache(self.click)
                echo_sys_exited(self.click)
                sys.exit(1)

            echo_local_fetched_cache_szie(self.click, len(results))

        results = self.dump_local._standardize_output(
            inputs, results, str(self.output_path), None, self.n_samples
        )
        self.generic_output_adapter._adapt_generic(
            json.dumps(results), str(self.output_path), self.model_id, DEFAULT_API_NAME
        )
        echo_redis_file_saved(self.click, str(self.output_path))
        sys.exit(1)
        return str(self.output_path)

    def _handle_local(self, inputs: list) -> str:
        self.dump_local.init_redis()
        missing_inputs = []
        if inputs:
            echo_redis_job_submitted(self.click, f"Input size: {len(inputs)}")
            results, _missing_inputs = self.dump_local.get_cached(
                self.model_id, inputs, self.dtype, cols=self.cols
            )
            echo_local_fetched_cache_szie(self.click, len(results))
            missing_inputs.extend(_missing_inputs)
            results = self.dump_local._standardize_output(
                inputs, results, str(self.output_path), None, self.n_samples
            )
            echo_redis_fetched_missed(self.click, len(results), len(missing_inputs))
            if len(results) + 1 >= len(inputs) and not missing_inputs:
                self.generic_output_adapter._adapt_generic(
                    json.dumps(results),
                    str(self.output_path),
                    self.model_id,
                    DEFAULT_API_NAME,
                )
                echo_redis_local_completed(self.click)
                echo_sys_exited(self.click)
                sys.exit(1)

            self.generic_output_adapter._adapt_generic(
                json.dumps(results),
                self.local_cache_csv_path,
                self.model_id,
                DEFAULT_API_NAME,
            )
        else:
            ns = self.n_samples if self.n_samples != -1 else "all"
            echo_redis_job_submitted(self.click, f"Sample size: {ns}")
            results, inputs = self.dump_local.fetch_all_cached(
                self.model_id, self.dtype
            )
            echo_local_fetched_cache_szie(self.click, len(results))
            if len(results) > self.n_samples:
                results = results[: self.n_samples]
            if len(results) < self.n_samples:
                self.n_samples = self.n_samples - len(results)
                results = self.dump_local._standardize_output(
                    inputs, results, str(self.output_path), None, self.n_samples
                )
                self.generic_output_adapter._adapt_generic(
                    json.dumps(results),
                    self.local_cache_csv_path,
                    self.model_id,
                    DEFAULT_API_NAME,
                )
        return results, missing_inputs

    def _get_none_size(self, res):
        size = 0
        for r in res:
            vals = list(r["output"].values())
            if "None" in vals[0]:
                size += 1
        return size

    def _submit_and_get_shards(self, inputs: list) -> list:
        if self.output_source == OutputSource.CACHE_ONLY:
            results, missing_input = self._handle_local(inputs)
            none_count = self._get_none_size(results)
            print(len(inputs), none_count)
            cache_size = abs(len(inputs) - none_count)
            echo_local_sample_warning(self.click, self.n_samples, cache_size)
            inputs = missing_input if len(missing_input) >= 1 else inputs

        s = self.n_samples if inputs is None else len(inputs)

        echo_small_sample_warning(self.click, s)
        self.request_id = str(uuid.uuid4())

        if (
            self.output_source in (OutputSource.CACHE_ONLY, OutputSource.CLOUD_ONLY)
            and inputs
        ):
            self.fetch_type = "filter"
            echo_uploading_inputs(self.click)
            pres = self.api.get_json(
                f"{INFERENCE_STORE_API_URL}/upload-destination",
                params={"modelid": self.model_id, "requestid": self.request_id},
            )
            tmp_path = self.files.create_temp_csv(inputs, self.input_adapter)
            self.files.upload_to_s3(pres, tmp_path)
            echo_upload_complete(self.click)
        echo_submitting_job(self.click, self.model_id)
        payload = {
            "requestid": self.request_id,
            "modelid": self.model_id,
            "fetchtype": self.fetch_type,
            "nsamples": self.n_samples,
            "dim": len(self.cols),
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

        return shards

    def _merge_shards(self, shards: list) -> str:
        self.files.merge_shards(
            shards,
            self.header,
            self.output_path,
            self.api,
            self.click,
            self.local_cache_csv_path,
        )
        echo_merged_saved(self.click, self.output_path)
        sys.exit(1)
        return str(self.output_path)
