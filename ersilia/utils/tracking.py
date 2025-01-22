import json
import os

from .exceptions_utils.throw_ersilia_exception import throw_ersilia_exception
from .session import get_session_dir, get_session_uuid

TRACKING_STUB = {
    "model_id": "",  # ID of the model
    "runs": 0,  # Number of runs
    "input_size": [],  # List of input sizes, where len(input_sizes) == runs
    "output_size": [],  # List of output sizes, where len(output_sizes) == runs
    "avg_input_size": [],  # List of average input sizes, where len(avg_input_sizes) == runs
    "avg_output_size": [],  # List of average output sizes, where len(avg_output_sizes) == runs
    "container_memory_perc": [],  # List of container memory percentages, where len(cpu_time_seconds) == runs
    "peak_container_memory_perc": [],  # List of peak container memory in MB, where len(peak_container_memory_MB) == runs
    "container_cpu_perc": [],  # List of container CPU utilisation percentages, where len(container_cpu_time_seconds) == runs
    "peak_container_cpu_perc": [],  # List of peak container CPU utilisation percentages, where len(peak_container_cpu_perc) == runs
    "nan_count_agg": [],  # List of NaN count, where len(nan_count_agg) == runs
    "mismatched_type_count": [],  # List of mismatched type count, where len(mismatched_type_agg) == runs
    "correct_shape": [],  # List of correct shape count, where len(correct_shape) == runs
    "warning_count": [],  # List of warning count, where len(warning_count) == runs
    "error_count": [],  # List of error count, where len(error_count) == runs
}

RUN_DATA_STUB = {  # This should not be used as is, always create a deep copy.
    "input_size": 0,
    "output_size": 0,
    "avg_input_size": 0,
    "avg_output_size": 0,
    "container_memory_perc": 0,
    "peak_container_memory_perc": 0,
    "container_cpu_perc": 0,
    "peak_container_cpu_perc": 0,
    "nan_count_agg": 0,
    "mismatched_type_count": 0,
    "correct_shape": False,
    "warning_count": 0,
    "error_count": 0,
}


@throw_ersilia_exception()
def init_tracking_summary(model_id):
    file_path = os.path.join(get_session_dir(), f"{get_session_uuid()}.json")
    try:
        with open(file_path, "w") as f:
            TRACKING_STUB["model_id"] = (
                model_id  # This won't affect models that run in different sessions, or if a new model is served in the same session.
            )
            json.dump(TRACKING_STUB, f)

    except FileNotFoundError:
        raise FileNotFoundError(
            "Could not create tracking summary file, check if session directory exists"
        )


@throw_ersilia_exception()
def update_tracking_summary(model_id, incoming_data):
    try:
        file_path = os.path.join(get_session_dir(), f"{get_session_uuid()}.json")
        with open(file_path, "r") as f:
            current = json.load(f)
        with open(file_path, "w") as f:
            current["runs"] += 1
            if current["model_id"] != model_id:
                raise ValueError("Model ID mismatch for given tracking data")
            for k, v in incoming_data.items():
                if k in current:
                    current[k].append(v)
                else:
                    raise KeyError(f"Key {k} not found in tracking data")
            json.dump(current, f)

    except FileNotFoundError:
        raise FileNotFoundError(
            "Could not update tracking summary file, check if session directory exists"
        )
