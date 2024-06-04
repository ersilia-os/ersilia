from datetime import datetime
import json
import csv
import statistics
from collections import defaultdict
import sys
import tracemalloc
import tempfile
import logging
import boto3
from botocore.exceptions import ClientError, NoCredentialsError
import os
import re
import requests

PERSISTENT_FILE_PATH = os.path.abspath("current_session.txt")
# Temporary path to log files until log files are fixed
TEMP_FILE_LOGS = os.path.abspath("")


def log_files_metrics(file):
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
        with open(file, "r") as file:
            line = None
            for line in file:
                if not re.match(r"^\d{2}.\d{2}.\d{2} \| ", line):
                    # continuation of log
                    if ersilia_error_flag:
                        # catch the error name if hinted by previous line
                        error_name = line.rstrip()
                        errors[error_name] += 1
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
                        errors["Unknown Ersilia exception class"] += 1
                        ersilia_error_flag = False
                    if misc_error_flag:
                        errors[error_name] += 1
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

        write_persistent_file(f"Error count: {error_count}")
        if len(errors) > 0:
            write_persistent_file(f"Breakdown by error types:")
            for error in errors:
                write_persistent_file(f"{error}: {errors[error]}")
        write_persistent_file(f"Warning count: {warning_count}")
    except (IsADirectoryError, FileNotFoundError):
        logging.warning("Unable to calculate metrics for log file: log file not found")
        

def open_persistent_file(model_id):
    """
    Opens a new persistent file, specifically for a run of model_id
    :param model_id: The currently running model
    """
    with open(PERSISTENT_FILE_PATH, "w") as f:
        f.write("Session started for model: {0}\n".format(model_id))
        
        
def write_persistent_file(contents):
    """
    Writes contents to the current persistent file. Only writes if the file actually exists.
    :param contents: The contents to write to the file.
    """

    # Only write to file if it already exists (we're meant to be tracking this run)
    if os.path.isfile(PERSISTENT_FILE_PATH):
        with open(PERSISTENT_FILE_PATH, "a") as f:
            f.write(f"{contents}\n")


def close_persistent_file():
    """
    Closes the persistent file, renaming it to a unique name.
    """

    # Make sure the file actually exists before we try renaming
    if os.path.isfile(PERSISTENT_FILE_PATH):
        log_files_metrics(TEMP_FILE_LOGS)

        new_file_path = os.path.join(
            os.path.dirname(PERSISTENT_FILE_PATH),
            datetime.now().strftime("%Y-%m-%d%_H-%M-%S.txt"),
        )
        os.rename(PERSISTENT_FILE_PATH, new_file_path)
        

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


def read_csv(file_path):
    """
    Reads a CSV file and returns the data as a list of dictionaries.

    :param file_path: Path to the CSV file.
    :return: A list of dictionaries containing the CSV data.
    """
    with open(file_path, mode='r') as file:
        reader = csv.DictReader(file)
        data = [row for row in reader]
    return data

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


class RunTracker:
    """
    This class will be responsible for tracking model runs. It calculates the desired metadata based on a model's
    inputs, outputs, and other run-specific features, before uploading them to AWS to be ingested
    to Ersilia's Splunk dashboard.
    """

    def __init__(self):
        self.time_start = None
        self.memory_usage_start = 0

    # function to be called before model is run
    def start_tracking(self):
        """
        Runs any code necessary for the beginning of the run.
        Currently necessary for tracking the runtime and memory usage of a run.
        """
        self.time_start = datetime.now()
        tracemalloc.start()
        self.memory_usage_start = tracemalloc.get_traced_memory()[0]

    def sample_df(self, df, num_rows, num_cols):
        """
        Returns a sample of the dataframe, with the specified number of rows and columns.
        """
        return df.sample(num_rows, axis=0).sample(num_cols, axis=1)

    def stats(self, result):
        """
        Stats function: calculates the basic statistics of the output file from a model. This includes the
        mode (if applicable), minimum, maximum, and standard deviation.

        :param result: The path to the model's output file.
        :return: A dictionary containing the stats for each column of the result.
        """

        data = read_csv(result)

        # drop first two columns (key, input)
        for row in data:
            row.pop('key', None)
            row.pop('input', None)
            
        # Convert data to a column-oriented format
        columns = defaultdict(list)
        for row in data:
            for key, value in row.items():
                columns[key].append(float(value))

        # Calculate statistics
        stats = {}
        for column, values in columns.items():
            column_stats = {}
            column_stats["mean"] = statistics.mean(values)
            try:
                column_stats["mode"] = statistics.mode(values)
            except statistics.StatisticsError:
                column_stats["mode"] = None
            column_stats["min"] = min(values)
            column_stats["max"] = max(values)
            column_stats["std"] = statistics.stdev(values) if len(values) > 1 else 0

            stats[column] = column_stats

        return stats
    
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
        This class is responsible for checking the types of the output file against the expected types.
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
            if not all(type_dict.get(dtype) == metadata["Output Type"][0] for dtype in types):
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
        :return: The peak memory usage in bytes.
        """

        # Compare memory between peak and amount when we started
        peak_memory = tracemalloc.get_traced_memory()[1] - self.memory_usage_start
        tracemalloc.stop()

        return peak_memory
    

    def track(self, input, result, meta):
        """
        Tracks the results after a model run.
        """
        json_dict = {}
        input_data = read_csv(input)
        result_data = read_csv(result)

        json_dict["input_data"] = input_data
        json_dict["result_data"] = result_data

        json_dict["meta"] = meta

        model_id = meta["metadata"].get("Identifier", "Unknown")
        json_dict["model_id"] = model_id

        time = datetime.now() - self.time_start
        json_dict["time_taken"] = str(time)

        # checking for mismatched types
        nan_count = get_nan_counts(result_data)
        json_dict["nan_count"] = nan_count

        json_dict["check_types"] = self.check_types(result_data, meta["metadata"])

        json_dict["stats"] = self.stats(result)

        json_dict["file_sizes"] = self.get_file_sizes(input_data, result_data)

        json_dict["peak_memory_use"] = self.get_peak_memory()

        # TODO: Call CDD Vault tracking and upload API success to splunk

        # log results to persistent tracking file
        json_object = json.dumps(json_dict, indent=4)
        write_persistent_file(json_object)

        # Upload run stats to s3
        upload_to_s3(json_dict)