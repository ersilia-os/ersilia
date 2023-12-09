from datetime import datetime
import json
import pandas as pd
import tracemalloc
import tempfile
import logging
import boto3
from botocore.exceptions import ClientError
import os
import re

PERSISTENT_FILE_PATH = os.path.abspath("current_session.txt")
# Temporary path to log files
TEMP_FILE_LOGS = os.path.abspath("")


def log_files_metrics(file):
    error_count = 0
    warning_count = 0

    ersilia_error_flag = False
    misc_error_flag = False
    error_name = ""
    errors = {}
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
                        error_name = line.split('| ERROR    | ')[1].rstrip()
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


def read_csv(file):
    # reads csv file and returns Pandas dataframe
    return pd.read_csv(file)


def read_json(result):
    data = json.load(result)
    return data


def open_persistent_file(model_id):
    with open(PERSISTENT_FILE_PATH, "w") as f:
        f.write("Session started for model: {0}\n".format(model_id))


def write_persistent_file(contents):
    # Only write to file if it already exists (we're meant to be tracking this run)
    if os.path.isfile(PERSISTENT_FILE_PATH):
        with open(PERSISTENT_FILE_PATH, "a") as f:
            f.write(f"{contents}\n")


def close_persistent_file():
    # Make sure the file actually exists before we try renaming
    if os.path.isfile(PERSISTENT_FILE_PATH):
        log_files_metrics(TEMP_FILE_LOGS)

        new_file_path = os.path.join(
            os.path.dirname(PERSISTENT_FILE_PATH),
            datetime.now().strftime("%Y-%m-%d%_H-%M-%S.txt"),
        )
        os.rename(PERSISTENT_FILE_PATH, new_file_path)


def upload_to_s3(json_dict, bucket="t4sg-ersilia", object_name=None):
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
        except ClientError as e:
            logging.error(e)
            return False
    return True


class RunTracker:
    """
    This class will be responsible for tracking model runs. It calculates the desired metadata based on a model's
    inputs, outputs, and other run-specific features, before uploading them to Ersilia's Splunk dashboard.

    NOTE: Currently, the Splunk connection is not set up. For now, we will print tracking results to the console.
    """

    def __init__(self):
        self.time_start = None
        self.memory_usage_start = 0

    # function to be called before model is run
    def start_tracking(self):
        self.time_start = datetime.now()
        tracemalloc.start()
        self.memory_usage_start = tracemalloc.get_traced_memory()[0]

    def sample_df(self, df, num_rows, num_cols):
        """
        Returns a sample of the dataframe, with the specified number of rows and columns.
        """
        return df.sample(num_rows, axis=0).sample(num_cols, axis=1)

    def stats(self, result):
        dat = read_csv(result)

        # drop first two columns (key, input)
        dat = dat.drop(["key", "input"], axis=1)

        # calculate and print statistics
        stats = {}
        for column in dat:
            column_stats = {}
            column_stats["mean"] = dat[column].mean()
            if len(dat[column].mode()) == 1:
                column_stats["mode"] = dat[column].mode().iloc[0]
            else:
                column_stats["mode"] = None
            column_stats["min"] = dat[column].min()
            column_stats["max"] = dat[column].max()
            column_stats["std"] = dat[column].std()

            stats[column] = column_stats

        return stats

    def get_file_sizes(self, input_df, output_df):
        input_size = input_df.memory_usage(deep=True).sum() / 1024
        output_size = output_df.memory_usage(deep=True).sum() / 1024

        input_avg_row_size = input_size / len(input_df)
        output_avg_row_size = output_size / len(output_df)

        return {
            "input_size": input_size,
            "output_size": output_size,
            "avg_input_size": input_avg_row_size,
            "avg_output_size": output_avg_row_size,
        }

    def check_types(self, resultDf, metadata):
        typeDict = {"float64": "Float", "int64": "Int"}
        count = 0

        # ignore key and input columns
        dtypesLst = resultDf.loc[:, ~resultDf.columns.isin(["key", "input"])].dtypes

        for i in dtypesLst:
            if typeDict[str(i)] != metadata["Output Type"][0]:
                count += 1

        if len(dtypesLst) > 1 and metadata["Output Shape"] != "List":
            print("Not right shape. Expected List but got Single")
            correct_shape = False
        elif len(dtypesLst) == 1 and metadata["Output Shape"] != "Single":
            print("Not right shape. Expected Single but got List")
            correct_shape = False
        else:
            print("Output is correct shape.")
            correct_shape = True

        print("Output has", count, "mismatched types.\n")

        return {"mismatched_types": count, "correct_shape": correct_shape}

    def get_peak_memory(self):
        # Compare memory between peak and amount when we started
        peak_memory = tracemalloc.get_traced_memory()[1] - self.memory_usage_start
        tracemalloc.stop()

        return peak_memory

    def track(self, input, result, meta):
        """
        Tracks the results after a model run.
        """
        json_dict = {}
        input_dataframe = read_csv(input)
        result_dataframe = read_csv(result)

        json_dict["input_dataframe"] = input_dataframe.to_dict()
        json_dict["result_dataframe"] = result_dataframe.to_dict()

        json_dict["meta"] = meta

        model_id = meta["metadata"].get("Identifier", "Unknown")
        json_dict["model_id"] = model_id

        time = datetime.now() - self.time_start
        json_dict["time_taken"] = str(time)

        # checking for mismatched types
        nan_count = result_dataframe.isna().sum()
        json_dict["nan_count"] = nan_count.to_dict()

        json_dict["check_types"] = self.check_types(result_dataframe, meta["metadata"])

        json_dict["stats"] = self.stats(result)

        json_dict["file_sizes"] = self.get_file_sizes(input_dataframe, result_dataframe)

        json_dict["peak_memory_use"] = self.get_peak_memory()

        # log results to persistent tracking file
        json_object = json.dumps(json_dict, indent=4)
        write_persistent_file(json_object)

        # Upload run stats to s3
        upload_to_s3(json_dict)
