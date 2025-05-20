import configparser
import csv
import json
import os
import re
import uuid
from datetime import datetime

import boto3
from botocore.exceptions import ClientError, NoCredentialsError

from ..io.output_logger import TabularResultLogger
from ..utils.exceptions_utils.throw_ersilia_exception import throw_ersilia_exception
from ..utils.exceptions_utils.tracking_exceptions import NoAWSCredentialsError
from ..utils.session import get_session_dir, get_session_uuid
from ..utils.system import SystemChecker
from .base import ErsiliaBase

TRACKING_BUCKET = "ersilia-models-runs"


class AwsConfig(ErsiliaBase):
    """
    This class is responsible for retrieving AWS credentials from the environment variables or the AWS config file.
    It checks for the presence of AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, and AWS_REGION in the environment variables.
    If not found, it looks for them in the AWS config file located at ~/.aws/credentials and ~/.aws/config.
    If the credentials are found, they are returned as a dictionary.
    """

    def __init__(self):
        ErsiliaBase.__init__(self, config_json=None, credentials_json=None)

    def get(self):
        """
        Get the AWS credentials from the environment variables or the AWS config file.

        Returns
        -------
        dict
            A dictionary containing the AWS credentials.
        """
        self.logger.debug(
            "Getting AWS credentials from environment variables or config and credentials files."
        )
        AWS_ACCESS_KEY_ID = os.environ.get("AWS_ACCESS_KEY_ID", None)
        AWS_SECRET_ACCESS_KEY = os.environ.get("AWS_SECRET_ACCESS_KEY", None)
        AWS_REGION = os.environ.get("AWS_REGION", None)
        if not AWS_ACCESS_KEY_ID:
            aws_config_path = os.path.expanduser("~/.aws/credentials")
            config = configparser.ConfigParser()
            config.read(aws_config_path)
            AWS_ACCESS_KEY_ID = config.get(
                "default", "aws_access_key_id", fallback=None
            ) or config.get("profile default", "aws_access_key_id", fallback=None)
        if not AWS_SECRET_ACCESS_KEY:
            aws_config_path = os.path.expanduser("~/.aws/credentials")
            config = configparser.ConfigParser()
            config.read(aws_config_path)
            AWS_SECRET_ACCESS_KEY = config.get(
                "default", "aws_secret_access_key", fallback=None
            ) or config.get("profile default", "aws_secret_access_key", fallback=None)
        if not AWS_REGION:
            aws_config_path = os.path.expanduser("~/.aws/config")
            config = configparser.ConfigParser()
            config.read(aws_config_path)
            AWS_REGION = config.get("default", "region", fallback=None) or config.get(
                "profile default", "region", fallback=None
            )
        data = {
            "aws_access_key_id": AWS_ACCESS_KEY_ID,
            "aws_secret_key_access": AWS_SECRET_ACCESS_KEY,
            "aws_region": AWS_REGION,
        }
        if data["aws_access_key_id"] is None:
            self.logger.warning(
                "AWS_ACCESS_KEY_ID not found in environment variables or config file."
            )
        else:
            self.logger.info("AWS_ACCESS_KEY_ID: ***")
        if data["aws_secret_key_access"] is None:
            self.logger.warning(
                "AWS_SECRET_ACCESS_KEY not found in environment variables or config file."
            )
        else:
            self.logger.info("AWS_SECRET_ACCESS_KEY: ***")
        if data["aws_region"] is None:
            self.logger.warning(
                "AWS_REGION not found in environment variables or config file."
            )
        else:
            self.logger.info("AWS_REGION: ***")
        return data

    def is_valid(self):
        """
        Validate access to AWS.

        Returns
        -------
        bool
            True if the configured access in the system is valid, False otherwise.
        """
        sts = boto3.client("sts")
        try:
            sts.get_caller_identity()
            self.logger.debug("AWS Account ID: ***")
            self.logger.debug("User ARN: ***")
            self.logger.debug("User ID: ***")
            return True
        except ClientError:
            self.logger.warning("Invalid AWS credentials.")
            return False


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

    def __init__(self, model_id, use_case, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.use_case = use_case
        self.aws_config = AwsConfig()
        self.time_start = None
        self.memory_usage_start = 0
        self.model_id = model_id

        self.session_folder = get_session_dir()
        self.logger.debug(
            "Run tracker folder of the session: {0}".format(self.session_folder)
        )

        self.lake_folder = os.path.join(self.session_folder, "lake")
        os.makedirs(self.lake_folder, exist_ok=True)
        self.logger.debug(
            "Run tracker folder of the lake: {0}".format(self.lake_folder)
        )

        self.tabular_result_logger = TabularResultLogger()

    @throw_ersilia_exception()
    def validate_aws_access(self):
        """
        Validate access to AWS.

        Returns
        -------
        bool
            True if AWS_ACCESS_KEY_ID and AWS_SECRET_KEY are available in the system
        """
        if self.aws_config.is_valid():
            self.logger.debug("AWS credentials are valid.")
        else:
            self.logger.warning("Invalid AWS credentials.")
            raise NoAWSCredentialsError()

    def log_files_metrics(self, file_log):
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
            res_dict["warning_count"] = warning_count
            return res_dict
        except (IsADirectoryError, FileNotFoundError):
            self.logger.warning(
                "Unable to calculate metrics for log file: log file not found"
            )

    def get_file_sizes(self, input_file, output_file):
        """
        Calculate the size of the input and output dataframes.

        Parameters
        ----------
        input_file : str
            File path containing the input data.
        output_file : str
            File path containing the output data.

        Returns
        -------
        dict
            Dictionary containing the input size, output size.
        """

        input_size = os.path.getsize(input_file) / (1024 * 1024)
        output_size = os.path.getsize(output_file) / (1024 * 1024)

        return {
            "input_size": input_size,
            "output_size": output_size,
        }

    def summarize_output(self, output_file):
        """
        This method summarizes the output of a model run
        Parameters
        ----------
        output_file : str
            The path to the output file.
        Returns
        -------
        data : dict
            A dictionary containing the summarized data.
        """
        data = self.tabular_result_logger.summary(output_file)
        return data

    def upload_to_s3(self, event_id):
        """
        Upload event information into an S3 bucket

        Parameters
        ----------
        event_id : list
            Event identifier.

        Returns
        -------
        bool
            True if uploading completed successfully, False otherwise.
        """
        bucket = TRACKING_BUCKET

        self.logger.debug("Connecting to ersilia-models-runs AWS S3 bucket")
        aws_config = self.aws_config.get()
        aws_access_key_id = aws_config["aws_access_key_id"]
        aws_secret_access_key = aws_config["aws_secret_key_access"]
        region_name = aws_config["aws_region"]

        s3_client = boto3.client(
            "s3",
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key,
            region_name=region_name,
        )

        self.logger.debug("Uploading event data to AWS S3 bucket: {0}".format(bucket))
        event_json = os.path.abspath(
            os.path.join(self.session_folder, "events", "json", f"{event_id}.json")
        )

        self.logger.debug("Event file JSON: {0}".format(event_json))

        with open(event_json, "r") as f:
            json_data = json.load(f)

        with open(event_json, "w") as f:
            json.dump(json_data, f)

        try:
            s3_client.upload_file(
                event_json,
                bucket,
                "events_json/{0}.json".format(event_id),
            )

        except NoCredentialsError:
            self.logger.error(
                "Unable to upload tracking data to AWS: Credentials not found"
            )
        except ClientError as e:
            self.logger.error(e)
            return False
        return True

    def get_country(self):
        """
        Get the country of the user.
        This method uses retrieves the country information.
        If the country information is not available, it defaults to 'Unknown'.

        Returns
        -------
        str
            The country of the user.
        """
        sc = SystemChecker()
        country, _ = sc.get_country()
        if country is None:
            self.logger.warning(
                "Unable to get country information. Defaulting to 'Unknown'."
            )
            country = "Unknown"
        self.logger.debug("Country: {0}".format(country))
        return country

    @throw_ersilia_exception()
    def create_event_data(self, event_id, input, output, metadata, time_seconds):
        """
        Track the results of a model run.

        This method collects relevant data for the run, updates the session file with the stats,
        and uploads the data to AWS if credentials are available.

        Parameters
        ----------
        input : str
            The input data used in the model run.
        output : str
            The output data in the form of a CSV file path.
        meta : dict
            The metadata of the model.
        container_metrics : dict
            The container metrics data.
        """
        run_data = {}

        size_info = self.get_file_sizes(input, output)

        current_log_file_path = os.path.join(get_session_dir(), "current.log")
        console_log_file_path = os.path.join(get_session_dir(), "console.log")
        error_and_warning_info_current_log = self.log_files_metrics(
            current_log_file_path
        )
        error_and_warning_info_console_log = self.log_files_metrics(
            console_log_file_path
        )
        run_data["input_size"] = (
            size_info["input_size"] if size_info["input_size"] else -1
        )
        run_data["output_size"] = (
            size_info["output_size"] if size_info["output_size"] else -1
        )
        run_data["error_count"] = (
            error_and_warning_info_current_log["error_count"]
            + error_and_warning_info_console_log["error_count"]
        )
        run_data["warning_count"] = (
            error_and_warning_info_current_log["warning_count"]
            + error_and_warning_info_console_log["warning_count"]
        )

        country = self.get_country()

        result_summary = self.summarize_output(output)

        data = {
            "session_id": get_session_uuid(),
            "event_id": event_id,
            "country": country,
            "date": datetime.now().strftime("%Y-%m-%d"),
            "use_case": self.use_case.capitalize(),
            "model_id": self.model_id,
            "slug": metadata["Slug"],
            "task": metadata["Task"],
            "subtask": metadata["Subtask"],
            "source_type": metadata["Source Type"],
            "output_consistency": metadata["Output Consistency"],
            "input_count": result_summary["num_inputs"],
            "output_dim": result_summary["output_dim"],
            "empty_cells_prop": result_summary["prop_empty_cells"],
            "empty_rows_count": result_summary["full_empty_rows"],
            "input_size_mb": run_data["input_size"],
            "output_size_mb": run_data["output_size"],
            "time_sec": time_seconds,
            "log_error_count": run_data["error_count"],
            "log_warning_count": run_data["warning_count"],
        }

        return data

    def track(self, input, output, metadata, time_seconds):
        """
        Track the model run and upload to S3 bucket.
        This method collects relevant data for the run, updates the session file with the stats,
        and uploads the data to AWS if credentials are available.
        Parameters
        ----------
        input : str
            The input data used in the model run.
        output : str
            The output data in the form of a CSV file path.
        metadata : dict
            The metadata of the model.
        Returns
        -------
        None
            This method does not return any value.
        """
        event_id = str(uuid.uuid4())
        data = self.create_event_data(
            event_id=event_id,
            input=input,
            output=output,
            metadata=metadata,
            time_seconds=time_seconds,
        )
        self.logger.debug("Event ID: {0}".format(event_id))
        events_folder = os.path.join(self.session_folder, "events")
        if not os.path.exists(events_folder):
            os.makedirs(events_folder)
        events_folder_json = os.path.join(events_folder, "json")
        if not os.path.exists(events_folder_json):
            os.makedirs(events_folder_json)
        event_file_json = os.path.join(events_folder_json, f"{event_id}.json")
        with open(event_file_json, "w") as f:
            json.dump(data, f, indent=4)
        self.logger.debug("Event data saved to {0}".format(event_file_json))
        events_folder_csv = os.path.join(events_folder, "csv")
        if not os.path.exists(events_folder_csv):
            os.makedirs(events_folder_csv)
        event_file_csv = os.path.join(events_folder_csv, f"{event_id}.csv")
        with open(event_file_csv, "w") as f:
            writer = csv.writer(f, delimiter=",")
            writer.writerow(data.keys())
            writer.writerow(data.values())
        self.logger.debug("Event data saved to {0}".format(event_file_csv))
        self.upload_to_s3(event_id)
