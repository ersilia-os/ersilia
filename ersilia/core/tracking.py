from datetime import datetime
import json
import pandas as pd
import tracemalloc
# from ersilia import cli
import logging
import boto3
from botocore.exceptions import ClientError
import os

class RunTracker:
    """
    This class will be responsible for tracking model runs. It calculates the desired metadata based on a model's
    inputs, outputs, and other run-specific features, before uploading them to Ersilia's Splunk dashboard.

    NOTE: Currently, the Splunk connection is not set up. For now, we will print tracking results to the console.
    """

    def __init__(self):
        self.time_start = None

    # function to be called before model is run
    def start_tracking(self):
        self.time_start = datetime.now()

    def read_csv(self, file):
        # reads csv file and returns Pandas dataframe
        return pd.read_csv(file)

    def track(self, input, result, meta):
        """
        Tracks the results after a model run.
        """

        print("Run input file:", input)
        print(self.read_csv(input))

        print("Run output file:", result)
        print(self.read_csv(result))

        print("Model metadata:", meta)

        model_id = meta["metadata"].get("Identifier", "Unknown")
        print("Model ID:", model_id)

        time = datetime.now() - self.time_start
        print("Time taken:", time)

    def log_to_console(self, data):
        print(f"\n{json.dumps(data)}\n")

    def read_json(self, result):
        data = json.load(result)
        self.log_to_console(result)
        return data
    
    def start(self):
        tracemalloc.start()
        self.time_start = tracemalloc.get_traced_memory()[0]
    
    def track_memory(self):
        peak_memory = tracemalloc.get_traced_memory()[1] - self.time_start
        print(f"Peak memory: {peak_memory}")
        tracemalloc.stop()


def write_file(dict):
    str = json.dump(dict)
    tmp = tempfile.NamedTemporaryFile()

    with open(tmp.name, 'w') as f:
        f.write(str)

def upload_file(file_name, bucket, object_name=None):
    """Upload a file to an S3 bucket

    :param file_name: File to upload
    :param bucket: Bucket to upload to
    :param object_name: S3 object name. If not specified then file_name is used
    :return: True if file was uploaded, else False
    """

    # If S3 object_name was not specified, use file_name
    if object_name is None:
        object_name = datetime.now().strftime("%Y-%m-%d-%H-%M-%S") + '-' + os.path.basename(file_name)

    # Upload the file
    s3_client = boto3.client('s3')
    try:
        response = s3_client.upload_file(file_name, bucket, object_name)
    except ClientError as e:
        logging.error(e)
        return False
    return True