import os
import re
import sys
import csv
import json
import time
import boto3
import shutil
import psutil
import logging
import requests
import tempfile
import types
import resource
import statistics
import tracemalloc
from .session import Session
from datetime import datetime
from datetime import timedelta
from .base import ErsiliaBase
from collections import defaultdict
from ..default import EOS, ERSILIA_RUNS_FOLDER
from ..utils.docker import SimpleDocker
from ..utils.csvfile import CsvDataLoader
from ..io.output_logger import TabularResultLogger
from botocore.exceptions import ClientError, NoCredentialsError


def flatten_dict(data):
    """
    This will flatten the nested dictionaries from the generator into a single-level dictionary,
    where keys from all levels are merged into one dictionary.

    :flat_dict: Result returned in a dictionary
    """
    flat_dict = {}
    for outer_key, inner_dict in data.items():
        for inner_key, value in inner_dict.items():
            flat_dict[inner_key] = value
    return flat_dict


def log_files_metrics(file_log, model_id):
    """
    This function will log the number of errors and warnings in the log files.

    :param file: The log file to be read
    :return: None (writes to file)
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
                    errors[error_name] += 1

        json_dict = {}
        json_dict["Error count"] = error_count

        if len(errors) > 0:
            json_dict["Breakdown by error types"] = {}
            for error in errors:
                json_dict["Breakdown by error types"][error] = errors[error]
        json_dict["Warning count"] = warning_count
        json_object = json.dumps(json_dict, indent=4)
        write_persistent_file(json_object, model_id)
    except (IsADirectoryError, FileNotFoundError):
        logging.warning("Unable to calculate metrics for log file: log file not found")


def get_persistent_file_path(model_id):
    """
    Construct the persistent file path.
    :param model_id: The currently running model
    :return: The path to the persistent file
    """
    return os.path.join(
        EOS, ERSILIA_RUNS_FOLDER, "session", model_id, "current_session.txt"
    )


def create_persistent_file(model_id):
    """
    Create the persistent path file
    :param model_id: The currently running model
    """

    persistent_file_dir = os.path.dirname(get_persistent_file_path(model_id))
    os.makedirs(persistent_file_dir, exist_ok=True)
    file_name = get_persistent_file_path(model_id)

    with open(file_name, "w") as f:
        f.write("Session started for model: {0}\n".format(model_id))


def check_file_exists(model_id):
    """
    Check if the persistent file exists.
    :param model_id: The currently running model
    :return: True if the file exists, False otherwise.
    """

    return os.path.isfile(get_persistent_file_path(model_id))


def write_persistent_file(contents, model_id):
    """
    Writes contents to the current persistent file. Only writes if the file actually exists.
    :param contents: The contents to write to the file.
    :param model_id: The currently running model
    """
    if check_file_exists(model_id):
        file_name = get_persistent_file_path(model_id)
        with open(file_name, "a") as f:
            f.write(f"{contents}\n")

    else:
        raise FileNotFoundError(
            f"The persistent file for model {model_id} does not exist. Cannot write contents."
        )


def close_persistent_file(model_id):
    """
    Closes the persistent file, renaming it to a unique name.
    :param model_id: The currently running model
    """
    if check_file_exists(model_id):
        file_name = get_persistent_file_path(model_id)
        file_log = os.path.join(EOS, "console.log")
        log_files_metrics(file_log, model_id)
        new_file_path = os.path.join(
            os.path.dirname(file_name),
            datetime.now().strftime("%Y-%m-%d_%H-%M-%S.txt"),
        )
        os.rename(file_name, new_file_path)

    else:
        raise FileNotFoundError(
            f"The persistent file for model {model_id} does not exist. Cannot close file."
        )


def upload_to_s3(json_dict, bucket="ersilia-tracking", object_name=None):
    """Upload a file to an S3 bucket

    :param json_dict: JSON object to upload
    :param bucket: Bucket to upload to
    :param object_name: S3 object name. If not specified then we generate a name based on the timestamp and model id.
    :return: True if file was uploaded, else False
    """

    # If S3 object_name was not specified, use file_name
    if object_name is None:
        object_name = (
            datetime.now().strftime("%Y-%m-%d-%H-%M-%S") + "-" + json_dict["model_id"]
        )

    # Dump JSON into a temporary file to upload
    json_str = json.dumps(json_dict, indent=4)
    tmp = tempfile.NamedTemporaryFile()

    with open(tmp.name, "w") as f:
        f.write(json_str)
        f.flush()

        # Upload the file
        s3_client = boto3.client("s3")
        try:
            s3_client.upload_file(tmp.name, bucket, f"{object_name}.json")
        except NoCredentialsError:
            logging.error(
                "Unable to upload tracking data to AWS: Credentials not found"
            )
        except ClientError as e:
            logging.error(e)
            return False
    return True


def upload_to_cddvault(output_df, api_key):
    """
    This function takes in the output dataframe from the model run and uploads the data to CDD vault.

    NOTE: Currently, this is simply a skeleton of what the final code should look like. The TODO details
    what the remaining changes should look like.

    :param output_df: The output dataframe from the model run
    :param api_key: The API key for CDD Vault's API
    :return: Whether the API call was successful
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
    Calculates the number of None values in each key of a list of dictionaries.

    :param data_list: List of dictionaries containing the data
    :return: Dictionary containing the count of None values for each key
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

    return nan_count


class RunTracker(ErsiliaBase):
    """
    This class will be responsible for tracking model runs. It calculates the desired metadata based on a model's
    inputs, outputs, and other run-specific features, before uploading them to AWS to be ingested
    to Ersilia's Splunk dashboard.
    """

    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.time_start = None
        self.memory_usage_start = 0
        self.model_id = model_id

        # Initialize folders
        self.ersilia_runs_folder = os.path.join(EOS, ERSILIA_RUNS_FOLDER)
        os.makedirs(self.ersilia_runs_folder, exist_ok=True)

        self.metadata_folder = os.path.join(self.ersilia_runs_folder, "metadata")
        os.makedirs(self.metadata_folder, exist_ok=True)

        self.lake_folder = os.path.join(self.ersilia_runs_folder, "lake")
        os.makedirs(self.lake_folder, exist_ok=True)

        self.logs_folder = os.path.join(self.ersilia_runs_folder, "logs")
        os.makedirs(self.logs_folder, exist_ok=True)

        self.tabular_result_logger = TabularResultLogger()

    #    TODO: see the following link for more details
    #    https://github.com/ersilia-os/ersilia/issues/1165?notification_referrer_id=NT_kwDOAsB0trQxMTEyNTc5MDIxNzo0NjE2NzIyMg#issuecomment-2178596998

    #    def stats(self, result):
    #        """
    #        Stats function: calculates the basic statistics of the output file from a model. This includes the
    #        mode (if applicable), minimum, maximum, and standard deviation.
    #        :param result: The path to the model's output file.
    #        :return: A dictionary containing the stats for each column of the result.
    #        """

    #        data = read_csv(result)

    # drop first two columns (key, input)
    #        for row in data:
    #            row.pop('key', None)
    #            row.pop('input', None)

    # Convert data to a column-oriented format
    #        columns = defaultdict(list)
    #        for row in data:
    #            for key, value in row.items():
    #                columns[key].append(float(value))

    # Calculate statistics
    #        stats = {}
    #        for column, values in columns.items():
    #            column_stats = {}
    #            column_stats["mean"] = statistics.mean(values)
    #            try:
    #                column_stats["mode"] = statistics.mode(values)
    #            except statistics.StatisticsError:
    #                column_stats["mode"] = None
    #            column_stats["min"] = min(values)
    #            column_stats["max"] = max(values)
    #            column_stats["std"] = statistics.stdev(values) if len(values) > 1 else 0
    #
    #            stats[column] = column_stats

    #        return stats

    def update_total_time(self, model_id, start_time):
        """
        Method to track and update the Total time taken by model.
        :Param model_id: The currently running model.
        :Param start_time: The start time of the running model.
        """

        end_time = time.time()
        duration = end_time - start_time
        if check_file_exists(model_id):
            file_name = get_persistent_file_path(model_id)
            with open(file_name, "r") as f:
                lines = f.readlines()

            updated_lines = []
            total_time_found = False

            for line in lines:
                if "Total time taken" in line and not total_time_found:
                    try:
                        total_time_str = line.split(":")[1].strip()
                        total_time = float(total_time_str)
                        total_time += duration
                        formatted_time = str(timedelta(seconds=total_time))
                        updated_lines.append(f"Total time taken: {formatted_time}\n")
                        total_time_found = True
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing 'Total time taken' value: {e}")
                else:
                    updated_lines.append(line)

            if not total_time_found:
                updated_lines.append(f"Total time taken: {formatted_duration}\n")

            new_content = "".join(updated_lines)
            with open(file_name, "w") as f:
                f.write(f"{new_content}\n")
        else:
            new_content = f"Total time: {formatted_duration}\n"
            with open(file_name, "w") as f:
                f.write(f"{new_content}\n")

    def get_file_sizes(self, input_file, output_file):
        """
        Calculates the size of the input and output dataframes, as well as the average size of each row.

        :input_file: Pandas dataframe containing the input data
        :output_file: Pandas dataframe containing the output data
        :return: dictionary containing the input size, output size, average input size, and average output size
        """

        input_size = sys.getsizeof(input_file) / 1024
        output_size = sys.getsizeof(output_file) / 1024

        input_avg_row_size = input_size / len(input_file)
        output_avg_row_size = output_size / len(output_file)

        return {
            "input_size": input_size,
            "output_size": output_size,
            "avg_input_size": input_avg_row_size,
            "avg_output_size": output_avg_row_size,
        }

    def check_types(self, result, metadata):
        """
        This method is responsible for checking the types of the output file against the expected types.
        This includes checking the shape of the output file (list vs single) and the types of each column.

        :param result: The output file
        :param metadata: The metadata dictionary
        :return: A dictionary containing the number of mismatched types and a boolean for whether the shape is correct
        """

        type_dict = {"float": "Float", "int": "Int"}
        count = 0

        # Collect data types for each column, ignoring "key" and "input" columns
        dtypes_list = {}
        for item in result:
            for key, value in item.items():
                if key not in ["key", "input"]:
                    if key not in dtypes_list:
                        dtypes_list[key] = set()
                    dtypes_list[key].add(type(value).__name__)

        mismatched_types = 0
        for column, types in dtypes_list.items():
            if not all(
                type_dict.get(dtype) == metadata["Output Type"][0] for dtype in types
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

        logging.info("Output has", count, "mismatched types.\n")

        return {"mismatched_types": count, "correct_shape": correct_shape}

    def get_peak_memory(self):
        """
        Calculates the peak memory usage of ersilia's Python instance during the run.
        :return: The peak memory usage in Megabytes.
        """

        usage = resource.getrusage(resource.RUSAGE_SELF)
        peak_memory_kb = usage.ru_maxrss
        peak_memory = peak_memory_kb / 1024
        return peak_memory

    def get_memory_info(self, process="ersilia"):
        """
        Retrieves the memory information of the current process
        """
        try:
            current_process = psutil.Process()
            process_name = current_process.name()
            cpu_times = current_process.cpu_times()

            if process_name != process:
                raise Exception(
                    f"Unexpected process. Expected: {process}, but got: {process_name}"
                )

            uss_mb = current_process.memory_full_info().uss / (1024 * 1024)
            total_cpu_time = sum(
                cpu_time
                for cpu_time in (
                    cpu_times.user,
                    cpu_times.system,
                    cpu_times.children_user,
                    cpu_times.children_system,
                    cpu_times.iowait,
                )
            )

            return uss_mb, total_cpu_time

        except psutil.NoSuchProcess:
            return "No such process found."
        except Exception as e:
            return str(e)

    def log_result(self, result):
        output_dir = os.path.join(self.lake_folder, self.model_id)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        file_name = os.path.join(output_dir, "{0}_lake.csv".format(self.model_id))
        tabular_result = self.tabular_result_logger.tabulate(result)
        if tabular_result is None:
            return
        with open(file_name, "w") as f:
            writer = csv.writer(f, delimiter=",")
            for r in tabular_result:
                writer.writerow(r)

    def log_meta(self, meta):
        output_dir = os.path.join(self.metadata_folder, self.model_id)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        file_name = os.path.join(output_dir, "{0}.json".format(self.model_id))
        with open(file_name, "w") as f:
            json.dump(meta, f)

    def log_logs(self):
        output_dir = os.path.join(self.logs_folder, self.model_id)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        file_name = os.path.join(output_dir, "{0}.log".format(self.model_id))
        session_file = os.path.join(EOS, "session.json")
        shutil.copyfile(session_file, file_name)

    def track(self, input, result, meta):
        """
        Tracks the results of a model run.
        """

        self.docker_client = SimpleDocker()
        self.data = CsvDataLoader()
        json_dict = {}

        if os.path.isfile(input):
            input_data = self.data.read(input)
        else:
            input_data = [{"SMILES": input}]

        # Create a temporary file to store the result if it is a generator
        if isinstance(result, types.GeneratorType):
            # Ensure EOS/tmp directory exists
            tmp_dir = os.path.join(EOS, "tmp")
            if not os.path.exists(tmp_dir):
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

        session = Session(config_json=self.config_json)
        model_id = meta["metadata"].get("Identifier", "Unknown")
        json_dict["model_id"] = model_id

        # checking for mismatched types
        nan_count = get_nan_counts(result_data)
        json_dict["nan_count"] = nan_count

        json_dict["check_types"] = self.check_types(result_data, meta["metadata"])

        json_dict["file_sizes"] = self.get_file_sizes(input_data, result_data)

        docker_info = (
            self.docker_client.container_memory(),
            self.docker_client.container_cpu(),
            self.docker_client.container_peak(),
        )

        json_dict["Docker Container"] = docker_info

        # Get the memory stats of the run processs
        peak_memory = self.get_peak_memory()
        total_memory, cpu_time = self.get_memory_info()

        # Update the session file with the stats
        session.update_peak_memory(peak_memory)
        session.update_total_memory(total_memory)
        session.update_cpu_time(cpu_time)
        self.log_logs()

        json_object = json.dumps(json_dict, indent=4)
        write_persistent_file(json_object, model_id)
        upload_to_s3(json_dict)

    # Log results, metadata, and session logs
    def log(self, result, meta):
        self.log_result(result)
        self.log_meta(meta)
        self.log_logs()
