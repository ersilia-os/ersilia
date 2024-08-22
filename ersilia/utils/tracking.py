import json
import os
from .session import get_session_dir, get_session_uuid
from .exceptions_utils.throw_ersilia_exception import throw_ersilia_exception



TRACKING_STUB = {
    "model_id": "", # ID of the model
    "runs": 0, # Number of runs
    "input_sizes": [], # List of input sizes, where len(input_sizes) == runs
    "output_sizes": [], # List of output sizes, where len(output_sizes) == runs
    "avg_input_sizes": [], # List of average input sizes, where len(avg_input_sizes) == runs
    "avg_output_sizes": [], # List of average output sizes, where len(avg_output_sizes) == runs
    "cpu_time_seconds": [], # List of CPU time in seconds, where len(cpu_time_seconds) == runs
    "container_memory_MB": [], # List of container memory in MB, where len(container_memory_MB) == runs
    "container_cpu_time_seconds": [], # List of container CPU time in seconds, where len(container_cpu_time_seconds) == runs
    "peak_container_memory_MB": [], # List of peak container memory in MB, where len(peak_container_memory_MB) == runs
    "nan_count_agg": [], # List of NaN count, where len(nan_count_agg) == runs
    "mismatched_type_agg": [], # List of mismatched type count, where len(mismatched_type_agg) == runs
    "correct_shape": [], # List of correct shape count, where len(correct_shape) == runs
    "warning_count": [], # List of warning count, where len(warning_count) == runs
    "error_count": [], # List of error count, where len(error_count) == runs
}

RUN_DATA_STUB = {  # This should not be used as is, always create a deep copy.
    "input_size": -1,
    "output_size": -1,
    "avg_input_size": -1,
    "avg_output_size": -1,
    "container_memory_perc": -1,
    "peak_container_memory_MB": -1,
    "container_cpu_perc": -1,
    "nan_count_agg": -1,
    "mismatched_type_count": -1,
    "correct_shape": False,
    "warning_count": -1,
    "error_count": -1
}

@throw_ersilia_exception
def init_tracking_summary(model_id):
    file_path = os.path.join(get_session_dir(), f"{get_session_uuid()}.json")
    try:
        with open(file_path, "w") as f:
            TRACKING_STUB["model_id"] = model_id # This won't affect models that run in different sessions, or if a new model is served in the same session.
            json.dump(TRACKING_STUB, f)
    
    except FileNotFoundError:
        raise FileNotFoundError("Could not create tracking summary file, check if session directory exists")


@throw_ersilia_exception
def update_tracking_summary(model_id, incoming_data):
    try:
        file_path = os.path.join(get_session_dir(), f"{get_session_uuid()}.json")
        with open(file_path, "rw") as f:
            current = json.load(f)
            current["runs"] += 1
            if current["model_id"] != model_id:
                raise ValueError("Model ID mismatch for given tracking data")
            for k,v in incoming_data:
                if k in current:
                    current[k].append(v)
                else:
                    raise KeyError(f"Key {k} not found in tracking data")

    except FileNotFoundError:
        raise FileNotFoundError("Could not update tracking summary file, check if session directory exists")