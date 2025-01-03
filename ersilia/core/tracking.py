import copy
import csv
import json
import os
import re
import resource
import sys
import tempfile
import types
from datetime import datetime

import boto3
import psutil
import requests
from botocore.exceptions import ClientError, NoCredentialsError
from loguru import logger as logging

from ..default import SESSION_JSON
from ..io.output_logger import TabularResultLogger
from ..utils.csvfile import CsvDataLoader
from ..utils.exceptions_utils.throw_ersilia_exception import throw_ersilia_exception
from ..utils.session import get_session_dir, get_session_uuid
from ..utils.tracking import (
    RUN_DATA_STUB,
    init_tracking_summary,
    update_tracking_summary,
)
from .base import ErsiliaBase
from .session import Session

AWS_ACCESS_KEY_ID = os.environ.get("AWS_ACCESS_KEY_ID")
AWS_SECRET_ACCESS_KEY = os.environ.get("AWS_SECRET_ACCESS_KEY")
AWS_REGION = os.environ.get("AWS_REGIOIN", "eu-central-1")
TRACKING_BUCKET = os.environ.get("TRACKING_BUCKET", "ersilia-models-runs")


def flatten_dict(data):
    """
    Flatten the nested dictionaries from the generator into a single-level dictionary.

    Parameters
    ----------
    data : dict
        The nested dictionary to flatten.

    Returns
    -------
    dict
        The flattened dictionary.
    """
    flat_dict = {}
    for outer_key, inner_dict in data.items():
        for inner_key, value in inner_dict.items():
            flat_dict[inner_key] = value
    return flat_dict


def log_files_metrics(file_log):
    """
    Log the number of errors and warnings in the log files.

    Parameters
    ----------
    file_log : str
        The log file to be read.

    Returns
    -------
    dict
        A dictionary containing the error count and warning count.
    """

    error_count = 0
    warning_count = 0

    ersilia_error_flag = False
    misc_error_flag = False
    error_name = ""
    errors = {}

    try:
        with open(file_log, "r") as file:
            line = None
            for line in file:
                if not re.match(r"^\d{2}.\d{2}.\d{2} \| ", line):
                    if ersilia_error_flag:
                        error_name = line.rstrip()
                        if error_name in errors:
                            errors[error_name] += 1
                        else:
                            errors[error_name] = 1
                        ersilia_error_flag = False
                        continue
                    elif misc_error_flag:
                        error_name += line.rstrip()
                        if len(error_name) > 100:
                            error_name = error_name[:97] + "..."
                            misc_error_flag = False
                else:
                    # encountering new logs
                    # make sure error flags are closed
                    if ersilia_error_flag:
                        errors["Unknown Ersilia exception class"] = (
                            errors.get("Unknown Ersilia exception class", 0) + 1
                        )
                        ersilia_error_flag = False
                    if misc_error_flag:
                        errors[error_name] = errors.get(error_name, 0) + 1
                        misc_error_flag = False
                    if "| ERROR" in line:
                        error_count += 1
                        # checking which type of errors
                        if "Ersilia exception class:" in line:
                            # combine this with the next line, usually EmptyOutputError or SourceCodeBaseInformationError
                            # the detailed message is long
                            ersilia_error_flag = True
                        else:
                            # other errors are pretty self-descriptive and short. Will cap by character
                            misc_error_flag = True
                            error_name = line.split("| ERROR    | ")[1].rstrip()
                    elif "| WARNING" in line:
                        warning_count += 1
            if line is not None:
                # in case last log is error
                # make sure error flags are closed
                if ersilia_error_flag:
                    errors["Unknown Ersilia exception class"] += 1
                if misc_error_flag:
                    errors[error_name] = 1 + errors.get(error_name, 0)

        res_dict = {}
        res_dict["error_count"] = error_count

        if len(errors) > 0:  # TODO We are not consuming this right now
            res_dict["error_details"] = {}
            for err in errors:
                res_dict["error_details"][err] = errors[err]
        res_dict["warning_count"] = warning_count
        return res_dict
    except (IsADirectoryError, FileNotFoundError):
        logging.warning("Unable to calculate metrics for log file: log file not found")


def serialize_session_json_to_csv(json_file, csv_file):
    """
    Serialize session JSON data to a CSV file.

    Parameters
    ----------
    json_file : str
        The path to the JSON file.
    csv_file : str
        The path to the CSV file.
    """
    with open(json_file, "r") as f:
        data = json.load(f)
        header = []
        values = []
        for k, v in data.items():
            header += [k]
            values += [v]
    with open(csv_file, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(values)


def serialize_tracking_json_to_csv(json_file, csv_file):
    """
    Serialize tracking JSON data to a CSV file.

    Parameters
    ----------
    json_file : str
        The path to the JSON file.
    csv_file : str
        The path to the CSV file.
    """
    with open(json_file, "r") as f:
        data = json.load(f)
        header = ["model_id"] + list(data.keys())[2:]  # Ignore model_id and runs
        num_rows = data["runs"]
        rows = []
        for i in range(num_rows):
            row = [data["model_id"]] + [data[k][i] for k in header[1:]]
            rows.append(row)
    with open(csv_file, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def upload_to_s3(model_id, metadata, bucket=TRACKING_BUCKET):
    """
    Upload a file to an S3 bucket.

    Parameters
    ----------
    model_id : str
        The identifier of the model.
    metadata : dict
        The metadata to upload.
    bucket : str, optional
        The S3 bucket to upload to, by default TRACKING_BUCKET.

    Returns
    -------
    bool
        True if the file was uploaded successfully, False otherwise.
    """

    s3_client = boto3.client(
        "s3",
        aws_access_key_id=AWS_ACCESS_KEY_ID,
        aws_secret_access_key=AWS_SECRET_ACCESS_KEY,
        region_name=AWS_REGION,
    )
    try:
        # Upload metadata to S3 by first writing to a temporary file
        tmp_metadata_file = tempfile.NamedTemporaryFile(mode="w", suffix=".json")
        with open(tmp_metadata_file.name, "w") as f:
            f.write(json.dumps(metadata, indent=4))

        s3_client.upload_file(
            tmp_metadata_file.name, bucket, f"metadata/{model_id}_metadata.json"
        )

        # Upload run output to S3
        sid = get_session_uuid()
        output_file_path = os.path.join(get_session_dir(), "lake", f"output_{sid}.csv")
        s3_client.upload_file(output_file_path, bucket, f"output/output_{sid}.csv")

        # Upload session info to S3
        session_json_path = os.path.join(get_session_dir(), SESSION_JSON)
        session_csv_path = session_json_path.split(".json")[0] + ".csv"
        serialize_session_json_to_csv(session_json_path, session_csv_path)
        s3_client.upload_file(session_csv_path, bucket, f"summary/session_{sid}.csv")
        os.remove(session_csv_path)

        # Upload tracking summary to S3
        tracking_json_path = os.path.join(
            get_session_dir(), f"{get_session_uuid()}.json"
        )
        s3_client.upload_file(tracking_json_path, bucket, f"tracking_raw/{sid}.json")
        tracking_csv_path = tracking_json_path.split(".json")[0] + ".csv"
        serialize_tracking_json_to_csv(tracking_json_path, tracking_csv_path)
        s3_client.upload_file(tracking_csv_path, bucket, f"tracking/{sid}.csv")
        os.remove(tracking_csv_path)

    except NoCredentialsError:
        logging.error("Unable to upload tracking data to AWS: Credentials not found")
    except ClientError as e:
        logging.error(e)
        return False
    return True


def upload_to_cddvault(output_df, api_key):
    """
    Upload the output dataframe from the model run to CDD Vault.

    Parameters
    ----------
    output_df : pd.DataFrame
        The output dataframe from the model run.
    api_key : str
        The API key for CDD Vault's API.

    Returns
    -------
    bool
        True if the API call was successful, False otherwise.
    """
    # We use the slurps API path to be able to bulk upload data
    url = "https://app.collaborativedrug.com/api/v1/vaults/<vault_id>/slurps"
    headers = {"CDD-Token": api_key}
    # TODO: Update project and header_mappings ids, as well as adding mappings for other
    # output columns if those are to be tracked as well.
    data = {
        "project": "",
        "autoreject": "true",
        "mapping_template": {
            "registration_type": "CHEMICAL_STRUCTURE",
            "header_mappings": [
                {
                    "header": {"name": "input", "position": 0},
                    "definition": {
                        "id": -1,
                        "type": "InternalFieldDefinition::MoleculeStructure",
                    },
                },
                {
                    "header": {"name": "time", "position": 1},
                    "definition": {
                        "id": -1,
                        "type": "InternalFieldDefinition::BatchFieldDefinition",
                    },
                },
            ],
        },
    }

    # Save output_df to a CSV of the correct format
    new_df = output_df[["input"]].copy()
    current_time = datetime.now().isoformat()

    new_df["time"] = current_time
    csv_file = tempfile.NamedTemporaryFile(mode="w", suffix=".csv")
    new_df.to_csv(csv_file.name, index=False)

    files = {"file": open(csv_file.name, "rb")}

    # Create and make API call
    response = requests.post(
        url, headers=headers, data={"json": json.dumps(data)}, files=files
    )

    if response.status_code == 200:
        return True
    else:
        logging.warning("API call to CDD Vault was Unsuccessful")
        return False


def get_nan_counts(data_list):
    """
    Calculate the number of NAN values in each key of a list of dictionaries.

    Parameters
    ----------
    data_list : list
        List of dictionaries containing the data.

    Returns
    -------
    int
        The count of NAN values for each key.
    """
    nan_count = {}

    # Collect all keys from data_list
    all_keys = set(key for item in data_list for key in item.keys())

    # Initialize nan_count with all keys
    for key in all_keys:
        nan_count[key] = 0

    # Count None values for each key
    for item in data_list:
        for key, value in item.items():
            if value is None:
                nan_count[key] += 1
    nan_count_agg = sum(nan_count.values())
    return nan_count_agg


class RunTracker(ErsiliaBase):
    """
    This class is responsible for tracking model runs. It calculates the desired metadata based on a model's
    inputs, outputs, and other run-specific features, before uploading them to AWS to be ingested
    to Ersilia's Splunk dashboard.

    Parameters
    ----------
    model_id : str
        The identifier of the model.
    config_json : dict
        Configuration in JSON format.
    """

    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.time_start = None
        self.memory_usage_start = 0
        self.model_id = model_id

        # Initialize folders
        self.session_folder = get_session_dir()

        self.lake_folder = os.path.join(self.session_folder, "lake")
        os.makedirs(self.lake_folder, exist_ok=True)

        self.tabular_result_logger = TabularResultLogger()

    def get_file_sizes(self, input_file, output_file):
        """
        Calculate the size of the input and output dataframes, as well as the average size of each row.

        Parameters
        ----------
        input_file : pd.DataFrame
            Pandas dataframe containing the input data.
        output_file : pd.DataFrame
            Pandas dataframe containing the output data.

        Returns
        -------
        dict
            Dictionary containing the input size, output size, average input size, and average output size.
        """

        input_size = sys.getsizeof(input_file) / 1024
        output_size = sys.getsizeof(output_file) / 1024

        try:
            input_avg_row_size = input_size / len(input_file)
            output_avg_row_size = output_size / len(output_file)
        except ZeroDivisionError:
            self.logger.warning(
                "Encountered a ZeroDivisionError. No data in input or output file"
            )
            input_avg_row_size = -1
            output_avg_row_size = -1

        return {
            "input_size": input_size,
            "output_size": output_size,
            "avg_input_size": input_avg_row_size,
            "avg_output_size": output_avg_row_size,
        }

    def check_types(self, result, metadata):
        """
        Check the types of the output file against the expected types.

        This method checks the shape of the output file (list vs single) and the types of each column.

        Parameters
        ----------
        result : list
            The output data.
        metadata : dict
            The metadata dictionary.

        Returns
        -------
        dict
            A dictionary containing the number of mismatched types and a boolean for whether the shape is correct.
        """

        type_dict = {"float": "Float", "int": "Int"}

        # Collect data types for each column, ignoring "key" and "input" columns
        dtypes_list = {}
        for item in result:
            for key, value in item.items():
                if key not in ["key", "input"]:
                    if key not in dtypes_list:
                        dtypes_list[key] = set()
                    dtypes_list[key].add(type(value).__name__)

        mismatched_types = 0
        for column, dtype_set in dtypes_list.items():
            if not all(
                type_dict.get(dtype) == metadata["Output Type"][0]
                for dtype in dtype_set
            ):
                mismatched_types += 1

        # Check if the shape is correct
        correct_shape = True
        if len(dtypes_list) > 1 and metadata["Output Shape"] != "List":
            logging.warning("Not right shape. Expected List but got Single")
            correct_shape = False
        elif len(dtypes_list) == 1 and metadata["Output Shape"] != "Single":
            logging.warning("Not right shape. Expected Single but got List")
            correct_shape = False

        logging.info(f"Output has {mismatched_types} mismatched types")

        return {"mismatched_types": mismatched_types, "correct_shape": correct_shape}

    def get_peak_memory(self):
        """
        Calculate the peak memory usage of Ersilia's Python instance during the run.

        Returns
        -------
        float
            The peak memory usage in Megabytes.
        """

        usage = resource.getrusage(resource.RUSAGE_SELF)
        peak_memory_kb = usage.ru_maxrss
        peak_memory = peak_memory_kb / 1024
        return peak_memory

    def get_memory_info(self):
        """
        Retrieve the memory information of the current process.

        Returns
        -------
        tuple
            A tuple containing the memory usage in MB and the total CPU time.
        """
        try:
            current_process = psutil.Process()
            cpu_times = current_process.cpu_times()

            uss_mb = current_process.memory_full_info().uss / (1024 * 1024)
            total_cpu_time = sum(
                cpu_time
                for cpu_time in (
                    cpu_times.user,
                    cpu_times.system,
                    cpu_times.children_user,
                    cpu_times.children_system,
                    # cpu_times.iowait,  # Is not platform agnostic
                )
            )

            return uss_mb, total_cpu_time

        except psutil.NoSuchProcess:
            logging.error("No such process found.")
            return "No such process found."
        except Exception as e:
            return str(e)

    def log_result(self, result):
        """
        Log the result of the model run.

        This method logs the result of the model run to a CSV file.

        Parameters
        ----------
        result : list
            The result data.
        """
        identifier = get_session_uuid()
        output_file = os.path.join(self.lake_folder, f"output_{identifier}.csv")
        tabular_result = self.tabular_result_logger.tabulate(
            result, identifier=identifier, model_id=self.model_id
        )
        if tabular_result is None:
            return
        with open(output_file, "a+") as f:
            writer = csv.writer(f, delimiter=",")
            for r in tabular_result:
                writer.writerow(r)

    @throw_ersilia_exception()
    def track(self, input, result, meta, container_metrics):
        """
        Track the results of a model run.

        This method collects relevant data for the run, updates the session file with the stats,
        and uploads the data to AWS if credentials are available.

        Parameters
        ----------
        input : str
            The input data used in the model run.
        result : str or Generator
            The output data in the form of a CSV file path or Generator from the model run.
        meta : dict
            The metadata of the model.
        container_metrics : dict
            The container metrics data.
        """
        # Set up requirements for tracking the run
        # self.docker_client = SimpleDocker()
        self.data = CsvDataLoader()
        run_data = copy.deepcopy(RUN_DATA_STUB)
        session = Session(config_json=self.config_json)
        model_id = meta["Identifier"]
        init_tracking_summary(model_id)

        if os.path.isfile(input):
            input_data = self.data.read(input)
        else:
            input_data = [{"input": input}]

        # Create a temporary file to store the result if it is a generator
        if isinstance(result, types.GeneratorType):
            tmp_dir = os.path.join(get_session_dir(), "tmp")
            os.makedirs(tmp_dir, exist_ok=True)

            # Create a temporary file to store the generator output
            temp_output_file = tempfile.NamedTemporaryFile(
                delete=False, suffix=".csv", dir=tmp_dir
            )

            flat_data_list = [flatten_dict(row) for row in result]
            if flat_data_list:
                header = list(flat_data_list[0].keys())
            temp_output_path = temp_output_file.name
            with open(temp_output_path, "w", newline="") as csvfile:
                csvWriter = csv.DictWriter(csvfile, fieldnames=header)
                csvWriter.writeheader()
                for flat_data in flat_data_list:
                    csvWriter.writerow(flat_data)
            result_data = self.data.read(temp_output_path)
            os.remove(temp_output_path)
        else:
            result_data = self.data.read(result)

        # Collect relevant data for the run
        nan_count = get_nan_counts(result_data)
        type_and_shape_info = self.check_types(result_data, meta)
        size_info = self.get_file_sizes(input_data, result_data)

        # peak_memory = self.docker_client.container_peak(self.model_id)
        current_log_file_path = os.path.join(get_session_dir(), "current.log")
        console_log_file_path = os.path.join(get_session_dir(), "console.log")
        error_and_warning_info_current_log = log_files_metrics(current_log_file_path)
        error_and_warning_info_console_log = log_files_metrics(console_log_file_path)
        run_data["input_size"] = (
            size_info["input_size"] if size_info["input_size"] else -1
        )
        run_data["output_size"] = (
            size_info["output_size"] if size_info["output_size"] else -1
        )
        run_data["avg_input_size"] = (
            size_info["avg_input_size"] if size_info["avg_input_size"] else -1
        )
        run_data["avg_output_size"] = (
            size_info["avg_output_size"] if size_info["avg_output_size"] else -1
        )
        run_data["container_cpu_perc"] = container_metrics["container_cpu_perc"]
        run_data["peak_container_cpu_perc"] = container_metrics["peak_cpu_perc"]
        run_data["container_memory_perc"] = container_metrics["container_memory_perc"]
        run_data["peak_container_memory_perc"] = container_metrics["peak_memory"]
        run_data["nan_count_agg"] = nan_count if nan_count else -1
        run_data["mismatched_type_count"] = type_and_shape_info["mismatched_types"]
        run_data["correct_shape"] = type_and_shape_info["correct_shape"]
        run_data["error_count"] = (
            error_and_warning_info_current_log["error_count"]
            + error_and_warning_info_console_log["error_count"]
        )
        run_data["warning_count"] = (
            error_and_warning_info_current_log["warning_count"]
            + error_and_warning_info_console_log["warning_count"]
        )

        update_tracking_summary(model_id, run_data)

        # Get the memory stats of the run processs
        peak_memory = self.get_peak_memory()
        total_memory, cpu_time = self.get_memory_info()

        # Update the session file with the stats
        session.update_peak_memory(peak_memory)
        session.update_total_memory(total_memory)
        session.update_cpu_time(cpu_time)
        self.log_result(result)

        if AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY:
            upload_to_s3(model_id=self.model_id, metadata=meta)
